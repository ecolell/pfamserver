# -*- coding: utf-8 -*-
from __future__ import print_function
from sqlalchemy import text
from application import app
from sqlalchemy import create_engine
from sqlalchemy_utils import database_exists, create_database
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
    "pfamA_reg_seed": "pfamA_acc,pfamseq_acc,seq_version,source,seq_start",
    "secondary_pfamseq_acc": "pfamseq_acc,secondary_acc",
    "pfamseq_disulphide": "pfamseq_acc,bond_start",
    "pfamseq_markup": "pfamseq_acc,auto_markup,residue",
    "pfamA_architecture": "pfamA_acc,auto_architecture",
    "interpro": "pfamA_acc",
    "pfamA_interactions": "pfamA_acc_A,pfamA_acc_B",
    "clan_membership": "pfamA_acc,clan_acc",
    "clan_alignment_and_relationship": "clan_acc",
    "clan_architecture": "clan_acc,auto_architecture",
    "pfamA_architecture": "pfamA_acc,auto_architecture",
    "dead_family": "pfamA_acc,pfamA_id",
    "dead_clan": "clan_acc,clan_id",
    "nested_locations": "pfamA_acc,nested_pfamA_acc,pfamseq_acc,seq_version,seq_start",
    "pdb_residue_data": "pdb_id,pfamseq_acc,pdb_seq_number,pfamseq_seq_number,pdb_res,pfamseq_res,pdb_insert_code,chain,serial,dssp_code",
    "pdb_image": "pdb_id",
    "proteome_regions": "pfamA_acc,ncbi_taxid",
    "complete_proteomes": "ncbi_taxid",
    "taxonomy": "ncbi_taxid",
    "pfamA2pfamA_scoop": "pfamA_acc_1,pfamA_acc_2",
    "pfamA_HMM": "pfamA_acc",
    "alignment_and_tree": "pfamA_acc,type",
}


tables = ['version', 'pfamA', 'pfamseq', 'uniprot']
tables += ['pfamA_reg_seed', 'uniprot_reg_full', 'pfamA_reg_full_significant',
           'pfamA_reg_full_insignificant']
tables += ['pfam_annseq', 'secondary_pfamseq_acc', 'evidence']
tables += ['markup_key', 'pfamseq_disulphide', 'other_reg', 'pfamseq_markup']
tables += ['architecture', 'pfamA_architecture']
tables += ['literature_reference', 'gene_ontology', 'pfamA_database_links',
           'interpro', 'pfamA_literature_reference', 'pfamA_interactions']
tables += ['clan_alignment_and_relationship', 'clan', 'clan_database_links',
           'clan_membership', 'clan_lit_ref', 'clan_architecture']
tables += ['dead_family', 'dead_clan']
tables += ['nested_locations']
tables += ['pdb', 'pdb_image', 'pdb_residue_data', 'pdb_pfamA_reg']
tables += ['proteome_regions', 'complete_proteomes', 'taxonomy',
           'ncbi_taxonomy']
tables += ['pfamA2pfamA_scoop', 'pfamA2pfamA_hhsearch']
tables += ['pfamA_HMM', 'alignment_and_tree']


class DatabaseConstructor(object):

    def __init__(self, version):
        self.version = version.version
        self.manager = version.manager
        self.backup_path = version.path + '/database_files/'
        self.url = app.config['SQLALCHEMY_DATABASE_URI'] + self.version
        if not database_exists(self.url):
            create_database(self.url)
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
                field_name = lambda r: any(map(lambda e: e in keys,
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


def init_db(version):
    constructor = DatabaseConstructor(version)
    constructor.construct()
