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
def db():
    """Database commands"""
    pass


tables = [
    'pfamA',
    'pfamseq',
    'uniprot',
    'pdb',
    'pdb_pfamA_reg',
    'uniprot_reg_full',
    'pfamA_reg_full_significant'
]


@db.group()
def structure():
    """Structure commands"""
    pass


def run(cmds, **kwargs):
    cmd = ' && '.join(['({:})'.format(c) for c in cmds]).format(**kwargs)
    click.echo("Running:\n{}".format(cmd))
    os.system(cmd)


def last_available_version():
    """Get the last available version."""
    url = 'http://ftp.ebi.ac.uk/pub/databases/Pfam/releases/?C=M;O=D'
    with closing(urlopen(url)) as conn:
        hrefs = [
            line for line in conn.readlines()
            if "href=\"Pfam" in line
        ]
        versions = [
            re.sub('^.+href="Pfam([0-9.]+).+\n$', r'\1', links)
            for links in hrefs
        ]
    return versions[0]


@structure.command()
@click.option('--version', '-v', 'version',
              type=click.STRING,
              multiple=False,
              default=last_available_version(),
              help='Version to install.')
@click.option('--ftp',
              is_flag=True,
              help='Force to use ftp.')
def build(version, ftp):
    """Get download commands."""
    protocol = 'ftp' if ftp else 'http'
    url = 'wget -c {protocol}://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam{version}/database_files/{file} -O ./Pfam{version}/{file}'
    files = [t + '.sql.gz' for t in tables]
    commands = [url.format(version=version, file=filename, protocol=protocol)
                for filename in files]
    commands = ['mkdir -p Pfam{version}'] + commands
    commands += [
        'gunzip -c Pfam{version}/{f_in} > Pfam{version}/{f_out}'.format(
            version=version,
            f_in=f,
            f_out=f[:-3]
        )
        for f in files]
    commands += [
        'echo "CREATE DATABASE IF NOT EXISTS {db_name}" | sudo mysql -u root',
        'cat Pfam{version}/*.sql | sudo mysql -u root {db_name}'
    ]
    db_name = 'Pfam' + version[:2] + '_' + version[-1:]
    run(commands, version=version, db_name=db_name)


@structure.command()
@click.option('--version', '-v', 'version',
              type=click.STRING,
              multiple=False,
              default=last_available_version(),
              help='Version to install.')
def clean(version):
    """Get download commands."""
    db_name = 'Pfam' + version[:2] + '_' + version[-1:]
    if click.confirm('Do you want to remove database {db_name}?'.format(db_name=db_name)):
        commands = [
            'echo "DROP DATABASE IF EXISTS {db_name}" | sudo mysql -u root',
        ]
        run(commands, version=version, db_name=db_name)
        click.echo('Ok done!')


@db.group()
def data():
    """Data commands"""
    pass


@data.command()
@click.option('--version', '-v', 'version',
              type=click.STRING,
              multiple=False,
              default=last_available_version(),
              help='Version to install.')
@click.option('--ftp',
              is_flag=True,
              help='Force to use ftp.')
def download(version, ftp):
    """Download files, continue on break (it could take more than a day)."""
    protocol = 'ftp' if ftp else 'http'
    url = 'wget -c {protocol}://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam{version}/database_files/{file} -O ./Pfam{version}/{file}'
    files = [t + '.txt.gz' for t in tables]
    commands = [
        'mkdir -p Pfam{version}'
    ]
    commands += [url.format(version=version, file=filename, protocol=protocol)
                 for filename in zip(files, tables)]
    click.echo('\n'.join(commands))
    run(commands, version=version)


@data.command()
@click.option('--version', '-v', 'version',
              type=click.STRING,
              multiple=False,
              default=last_available_version(),
              help='Version to install.')
def load(version):
    """Load the data into the database (be aware to apply into a clean database)."""
    files = [(t + '.txt.gz') for t in tables]
    root = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '../..'))
    click.echo(root)
    commands = [
        'mkdir -p Pfam{version}'
    ]
    loader = ['sudo mysql -u root {db_name} -e "LOAD DATA LOCAL INFILE \'{root}/Pfam{version}/{f_out}\'',
              'INTO TABLE {table} CHARACTER SET latin1 COLUMNS TERMINATED BY \'\\t\' LINES TERMINATED BY \'\\n\';"']
    command = ' && '.join(
        [
            'gunzip -c Pfam{version}/{f_in} > Pfam{version}/{f_out}',
            ' '.join(loader),
            'rm Pfam{version}/{f_out}'
        ]
    )
    db_name = 'Pfam' + version[:2] + '_' + version[-1:]
    commands += [command.format(version=version, f_in=filename, f_out=filename[:-3], root=root, table=table, db_name=db_name)
                 for (filename, table) in zip(files, tables)]
    click.echo('\n'.join(commands))
    run(commands, version=version)
