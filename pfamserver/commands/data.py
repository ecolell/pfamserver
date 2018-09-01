# -*- coding: utf-8 -*-
from __future__ import absolute_import
from __future__ import print_function
from __future__ import unicode_literals

import click
click.disable_unicode_literals_warning = True

from urllib import urlopen
from contextlib import closing
import re
import os


@click.group()
def data():
    """Data commands"""
    pass


def run(cmds, **kwargs):
    cmd = ' && '.join(['({:})'.format(c) for c in cmds]).format(**kwargs)
    click.echo("Running:\n{}".format(cmd))
    os.system(cmd)


def last_available_version():
    """Get the last available version."""
    url = 'http://ftp.ebi.ac.uk/pub/databases/Pfam/releases/?C=M;O=D'
    with closing(urlopen(url)) as conn:
        hrefs = filter(lambda l: "href=\"Pfam" in l, conn.readlines())
        versions = map(lambda l: re.sub('^.+href="Pfam([0-9.]+).+\n$',
            r'\1', l), hrefs)
    return versions[0]


@data.command()
@click.option('--version', '-v', 'version',
              type=click.STRING,
              multiple=False,
              default=last_available_version(),
              help='Version to install.')
@click.option('--ftp',
              is_flag=True,
              help='Force to use ftp.')
def download_commands(version, ftp):
    """Get download commands."""
    protocol = 'ftp' if ftp else 'http'
    url = 'wget -c {protocol}://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam{version}/database_files/{file} -O ./Pfam{version}/{file}'
    files = [
        'pdb_pfamA_reg.sql.gz',
        'pdb_pfamA_reg.txt.gz',
        'pdb.sql.gz',
        'pdb.txt.gz',
        'pfamA_reg_full_significant.sql.gz',
        'pfamA_reg_full_significant.txt.gz',
        'pfamA.sql.gz',
        'pfamA.txt.gz',
        'pfamseq.sql.gz',
        'pfamseq.txt.gz',
        'uniprot_reg_full.sql.gz',
        'uniprot_reg_full.txt.gz',
        'uniprot.sql.gz',
        'uniprot.txt.gz'
    ]
    commands = [url.format(version=version, file=filename, protocol=protocol)
                for filename in files]
    click.echo('\n'.join(commands))
    run([
        'mkdir -p Pfam{version}'
    ], version=version)

# verify old process to create the required primary keys and undefined indexes to speed up the queries.
# sudo mysql -u root -p pfamserver -e "LOAD DATA LOCAL INFILE '/home/eloy/version/git/pfamserver/dumps/pfamA.txt'
# INTO TABLE pfamA CHARACTER SET latin1 COLUMNS TERMINATED BY '\t' LINES TERMINATED BY '\n';"
