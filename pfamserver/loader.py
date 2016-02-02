from sqlalchemy import text
from application import app
from database import db
import sqlparse
import os
import gzip
import shutil
import re


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


patched_tables = {
    "version": "pfam_release",
    "pfamA_reg_seed": "pfamA_acc",
    "pfam_annseq": "pfamseq_acc",
    "secondary_pfamseq_acc": "pfamseq_acc",
    "pfamseq_disulphide": "pfamseq_acc",
    "pfamseq_markup": "pfamseq_acc,auto_markup",
    "pfamA_architecture": "pfamA_acc,auto_architecture",
    "gene_ontology": "pfamA_acc,go_id",
    "pfamA_database_links": "pfamA_acc,db_id",
    "interpro": "pfamA_acc,interpro_id",
    "pfamA_interactions": "pfamA_acc_A,pfamA_acc_B",
    "pfamA_literature_reference": "pfamA_acc,auto_lit",
    "clan_membership": "pfamA_acc,clan_acc",
    "clan_alignment_and_relationship": "clan_clan_acc",
    "clan_database_links": "clan_acc,db_id",
    "clan_lit_ref": "clan_acc,auto_lit",
    "clan_architecture": "clan_acc,auto_architecture",
    "pfamA_architecture": "pfamA_acc,auto_architecture",
    "dead_family": "pfamA_acc,pfamA_id",
    "dead_clan": "clan_acc, clan_id",
    "nested_locations": "pfamA_acc,nested_pfamA_acc,pfamseq_acc",
    "pdb_residue_data": "pdb_id,pfamseq_acc",
    "pdb_image": "pdb_id,pdb_pdb_id",
    "proteome_regions": "pfamA_acc,ncbi_taxid",
    "complete_proteomes": "ncbi_taxid",
    "taxonomy": "ncbi_taxid",
    "pfamA2pfamA_scoop": "pfamA_acc_1,pfamA_acc_2",
    "pfamA2pfamA_hhsearch": "pfamA_acc_1,pfamA_acc_2",
    "pfamA_HMM": "pfamA_acc",
    "alignment_and_tree": "pfamA_acc",
}


@decompress
def load_sql(table_name):
    with app.open_resource('{:}{:}.sql'.format(backup_path, table_name),
                           mode='r') as f:
        # TODO: If not primary_key it should use a dictionary to modify the
        # files.
        cmds = f.read().decode('unicode_escape').encode('ascii', 'ignore')
        cmds = cmds.replace('NOT NULL', '')
        cmds = '\n'.join(filter(lambda c: c and c[0] not in ['/', '-'],
                              cmds.split('\n')))
        if 'PRIMARY KEY' not in cmds:
            if table_name in patched_tables:
                cmds = cmds.replace(') ENGINE',
                                    ', PRIMARY KEY ({:})) ENGINE'.format(
                                        patched_tables[table_name]))
        cmds = sqlparse.split(cmds)
        #map(execute, cmds)


@decompress
def load_txt(table_name):
    backup_filename = '{:}{:}.txt'.format(backup_path, table_name)

    cmd = ("LOAD DATA INFILE '{:}' INTO TABLE {:} COLUMNS TERMINATED BY '\t' "
           "LINES TERMINATED BY '\n'")
    cmd = cmd.format(backup_filename, table_name)
    #execute(cmd)


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


def init_db():
    map(load_sql, tables)
    #map(load_txt, tables)
