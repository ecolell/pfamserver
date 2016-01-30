from sqlalchemy import text
from config import app, db
import sqlparse
import os
import gzip
import shutil


backup_path = os.path.dirname(os.path.abspath(__file__))
backup_path += '/database_files/'


def decompress(function):
    extension = function.__name__[len('load_'):]

    def wrapper(table_name):
        destiny = backup_path + '{:}.{:}'.format(table_name, extension)
        source = destiny + '.gz'
        with gzip.open(source, 'rb') as f_in:
            with open(destiny, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        result = function(table_name)
        os.remove(destiny)
        return result
    return wrapper


def execute(command):
    try:
        db.engine.execute(text(command))
    except Exception, e:
        print e


@decompress
def load_sql(table_name):
    with app.open_resource('{:}{:}.sql'.format(backup_path, table_name),
                           mode='r') as f:
        cmds = f.read().decode('unicode_escape').encode('ascii', 'ignore')
        cmds = cmds.replace('NOT NULL', '')
        cmds = map(lambda c: c.split(';')[0], sqlparse.split(cmds))
        map(execute, cmds)


@decompress
def load_txt(table_name):
    backup_filename = '{:}{:}.txt'.format(backup_path, table_name)

    cmd = ("LOAD DATA INFILE '{:}' INTO TABLE {:} COLUMNS TERMINATED BY '\t' "
           "LINES TERMINATED BY '\n'")
    cmd = cmd.format(backup_filename, table_name)
    execute(cmd)


def init_db():
    tables = ['version', 'pfamA', 'pfamseq', 'uniprot']
    map(load_sql, tables)
    map(load_txt, tables)
