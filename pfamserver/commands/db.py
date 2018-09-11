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
from pfamserver.commands.unused_columns import unused_columns


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
    """Clean database."""
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
              help='Version to download.')
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
                 for filename in files]
    click.echo('\n'.join(commands))
    run(commands, version=version)


@data.command()
@click.option('--version', '-v', 'version',
              type=click.STRING,
              multiple=False,
              default=last_available_version(),
              help='Version to load.')
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


@data.command()
@click.option('--version', '-v', 'version',
              type=click.STRING,
              multiple=False,
              default=last_available_version(),
              help='Version to measure size.')
def size(version):
    """Get database size."""
    db_name = 'Pfam' + version[:2] + '_' + version[-1:]
    query = 'SELECT table_schema, ' \
            'ROUND(SUM(data_length + index_length) / 1024 / 1024 / 1024, 1) \'DB Size in GB\' ' \
            'FROM information_schema.tables ' \
            'GROUP BY table_schema HAVING table_schema=\'{db_name}\';'.format(db_name=db_name)
    commands = ['sudo mysql -u root {db_name} -e "{query}"'.format(query=query, db_name=db_name)]
    run(commands)


@data.command()
@click.option('--version', '-v', 'version',
              type=click.STRING,
              multiple=False,
              default=last_available_version(),
              help='Version to shrink.')
def shrink(version):
    """Shrink the data into the database removing unused columns."""
    db_name = 'Pfam' + version[:2] + '_' + version[-1:]
    query = 'SELECT \'{table}->{column}\'; ' \
            'set @exist_Check := ( ' \
            '   select count(*) from information_schema.columns ' \
            '   where table_name=\'{table}\' ' \
            '   and column_name=\'{column}\' ' \
            '   and table_schema=database()  ' \
            ') ; ' \
            'set @sqlstmt := if(@exist_Check>0,\'ALTER TABLE {table} DROP COLUMN {column};\' , \'select 1\') ; '\
            'prepare stmt from @sqlstmt ; '\
            'execute stmt ;'
    db_name = 'Pfam' + version[:2] + '_' + version[-1:]
    superquery = [
        query.format(db_name=db_name, table=table, column=column)
        for (table, columns) in unused_columns.items() for column in columns
    ]
    queries = [
        'UPDATE uniprot SET created=updated;',
        'UPDATE pfamseq SET created=updated;',
        'UPDATE pfamA SET created=updated;',
        ''.join(superquery)
    ]
    commands = [
        'sudo mysql -u root {db_name} -e "{query}"'.format(
            query=q,
            db_name=db_name
        )
        for q in queries
    ]
    run(commands)


@db.command()
@click.option('--version', '-v', 'version',
              type=click.STRING,
              multiple=False,
              default=last_available_version(),
              help='Version to dump.')
def dump(version):
    """Dump the database into a sql file."""
    db_name = 'Pfam' + version[:2] + '_' + version[-1:]
    commands = [
        'mkdir -p data',
        'sudo mysqldump --databases {db_name} > ./data/pfam{version}.sql'
    ]
    run(commands, db_name=db_name, version=version)


@db.command()
@click.option('--version', '-v', 'version',
              type=click.STRING,
              multiple=False,
              default=last_available_version(),
              help='Version to dump.')
def pack_dump(version):
    """Pack the dump of the database into a bzip2 file."""
    commands = [
        'bzip2 -c ./data/pfam{version}.sql > ./data/pfam{version}.sql.bz2'
    ]
    run(commands, version=version)


@db.group()
def shrinked():
    """Shrinked data commands"""
    pass


versions = {
    '31.0': '1S0jrbULH9ZlXrDsxQPYms29REuYDNjDK'
}


@shrinked.command()
@click.option('--version', '-v', 'version',
              type=click.STRING,
              multiple=False,
              default=last_available_version(),
              help='Version to install.')
def download(version):
    """Download shrinked file, continue on break (it could take a few hours)."""
    if version not in versions:
        click.echo('Version ' + version + ' is not available as a skrinked DB.')
        return None
    filename = './data/pfam' + version + '.sql.bz2'
    commands = [
        'mkdir -p data',
        'bash ./pfamserver/commands/shrinked_downloader.sh {id} {filename}'
    ]
    click.echo('\n'.join(commands))
    run(commands,
        id=versions[version],
        filename=filename)


@shrinked.command()
@click.option('--version', '-v', 'version',
              type=click.STRING,
              multiple=False,
              default=last_available_version(),
              help='Version to install.')
def install(version):
    """Install shrinked file, continue on break (it could take a few hours)."""
    filename = './data/pfam' + version + '.sql.bz2'
    data_filename = './data/pfam' + version + '.sql'
    if not os.path.isfile(filename):
        click.echo('Version ' + version + ' wasn\'t downloaded as a shrinked DB.')
        return None
    commands = [
        'bzip2 -d {filename} > {data_filename}',
        'cat {data_filename} | sudo mysql -u root'
    ]
    click.echo('\n'.join(commands))
    run(commands,
        filename=filename,
        data_filename=data_filename)
