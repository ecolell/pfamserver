# -*- coding: utf-8 -*-
from __future__ import print_function
from sqlalchemy import text
from application import app
from sqlalchemy import create_engine
import sqlparse
import os
import gzip
import shutil
from core import Manager
import re
import warnings
import json


def decompress(function):
    extension = function.__name__[len('load_'):]

    def gz_wrapper(*args, **kwargs):
        constructor = args[0]
        table_name = args[1]
        destiny = constructor.backup_path + '{:}.{:}'.format(table_name, extension)
        source = destiny + '.gz'
        with gzip.open(source, 'rb') as f_in:
            with open(destiny, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        result = function(*args, **kwargs)
        os.remove(destiny)
        return result
    return gz_wrapper, extension


patched_pk = {
    "version": "number_families",
}


tables = ['version', 'pfamA', 'pfamseq', 'uniprot']
tables += ['uniprot_reg_full', 'pfamA_reg_full_significant']
tables += ['pdb', 'pdb_pfamA_reg']


class DatabaseConstructor(object):

    def __init__(self, version):
        self.version = version.version
        self.manager = version.manager
        self.backup_path = version.path + '/database_files/'
        self.url = app.config['SQLALCHEMY_DATABASE_URI'] + self.version
        self.manager.prepare_database(self.url)
        self.engine = create_engine(self.url,
                                    pool_size=20,
                                    max_overflow=100)

    def construct(self):
        map(self.load_sql, tables)
        map(self.load_txt, tables)

    def execute(self, command):
        self.engine.execute(text(command))

    @Manager.milestone
    @decompress
    def load_sql(self, table_name):
        filename = '{:}{:}.sql'.format(self.backup_path, table_name)
        print("\t\tCreating {:} structure into the database.".format(table_name))
        with app.open_resource(filename, mode='r') as f:
            cmds = f.read().decode('unicode_escape').encode('ascii', 'ignore')
            cmds = cmds.replace('NOT NULL', '').split('\n')
            keys = []
            if table_name in patched_pk:
                keys = patched_pk[table_name].split(',')

                def field_name(r):
                    return any(map(lambda e: e in keys,
                                   re.findall('`[^`]*`', r)))
                cmds = map(lambda r: r if field_name(r) else
                           r.replace('DEFAULT NULL', ''),
                           cmds)
            cmds = filter(lambda c: c and c[0] not in ['/', '-'],
                          cmds)
            cmds = ''.join(cmds)
            cmds = cmds.replace("`created` datetime DEFAULT NULL",
                                "`created` datetime NULL DEFAULT NULL")
            if 'PRIMARY KEY' not in cmds:
                if table_name in patched_pk:
                    keys = map(lambda k: "`" + k + "`", keys)
                    cmds = cmds.replace(') ENGINE',
                                        ', PRIMARY KEY ({:})) ENGINE'.format(
                                            ','.join(keys)))
            cmds = sqlparse.split(cmds)
            with warnings.catch_warnings():
                warnings.filterwarnings('ignore', 'unknown table')
                warnings.filterwarnings('ignore', 'Duplicate index',
                                        append=True)
                map(self.execute, cmds)
            return True

    @Manager.milestone
    @decompress
    def load_txt(self, table_name):
        print("\t\tLoading {:} data to the database.".format(table_name))
        backup_filename = '{:}{:}.txt'.format(self.backup_path, table_name)

        cmd = ("LOAD DATA INFILE '{:}' INTO TABLE {:} CHARACTER SET latin1 "
               "COLUMNS TERMINATED BY '\t' LINES TERMINATED BY '\n'")
        cmd = cmd.format(backup_filename, table_name)
        self.execute("SET SESSION sql_mode='ALLOW_INVALID_DATES'")
        with warnings.catch_warnings():
            warnings.filterwarnings('ignore', 'Incorrect integer')
            warnings.filterwarnings('ignore', 'Incorrect decimal', append=True)
            warnings.filterwarnings('ignore', 'Data truncated', append=True)
            self.execute(cmd)
        return True

