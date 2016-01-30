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
    tables =  ['version', 'pfamA', 'pfamseq', 'uniprot']
    tables += ['pfamA_reg_seed', 'uniprot_reg_full', 'pfamA_reg_full_significant', 'pfamA_reg_full_insignificant']
    tables += ['pfam_annseq', 'secondary_pfamseq_acc', 'evidence']
    tables += ['markup_key', 'pfamseq_disulphide', 'other_reg', 'pfamseq_markup']
    tables += ['architecture', 'pfamA_architecture']
    tables += ['literature_reference', 'gene_ontology', 'pfamA_database_links', 'interpro', 'pfamA_literature_reference', 'pfamA_interactions']
    tables += ['clan_alignment_and_relationship', 'clan', 'clan_database_links', 'clan_membership', 'clan_lit_ref', 'clan_architecture']
    tables += ['dead_family', 'dead_clan']
    tables += ['nested_locations']
    tables += ['pdb', 'pdb_image', 'pdb_residue_data', 'pdb_pfamA_reg']
    tables += ['proteome_regions', 'complete_proteomes', 'taxonomy', 'ncbi_taxonomy']
    tables += ['pfamA2pfamA_scoop', 'pfamA2pfamA_hhsearch']
    tables += ['pfamA_HMM', 'alignment_and_tree']
    map(load_sql, tables)
    map(load_txt, tables)
