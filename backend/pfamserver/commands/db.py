import os

import click

click.disable_unicode_literals_warning = True


@click.group()
def db():
    """Database commands"""
    pass


tables = [
    "pfamA",
    "pfamseq",
    "uniprot",
    "pdb",
    "pdb_pfamA_reg",
    "uniprot_reg_full",
    "pfamA_reg_full_significant",
]


@db.group()
def data():
    """Data commands"""
    pass


@data.command()
@click.option(
    "--version",
    "-v",
    "version",
    type=click.STRING,
    multiple=False,
    default="35.0",
    help="Version to build the cache.",
)
def build_cache(version):
    """Build cache tables to improve the performance of some queries."""
    db_name = "Pfam" + version[:2] + "_" + version[-1:]
    commands = [
        'cat /home/pfamserver/stage/pfamserver/commands/db_build_cache.sql | mysql -u root -h db --password="root" {db_name}'.format(
            db_name=db_name
        )
    ]
    # run(commands)


@db.command()
@click.option(
    "--version",
    "-v",
    "version",
    type=click.STRING,
    multiple=False,
    default="35.0",
    help="Version to dump.",
)
def dump(version):
    """Dump the database into a sql file."""
    db_name = "Pfam" + version[:2] + "_" + version[-1:]
    commands = [
        "mkdir -p mysql",
        "sudo mysqldump --databasecat  | s {db_name} > ./mysql/pfam{version}.sql",
    ]
    # run(commands, db_name=db_name, version=version)


@db.command()
@click.option(
    "--version",
    "-v",
    "version",
    type=click.STRING,
    multiple=False,
    default="35.0",
    help="Version to dump.",
)
def pack_dump(version):
    """Pack the dump of the database into a bzip2 file."""
    commands = ["bzip2 -c ./mysql/pfam{version}.sql > ./mysql/pfam{version}.sql.bz2"]
    run(commands, version=version)


@db.group()
def shrinked():
    """Shrinked mysql commands"""
    pass


versions = {
    "31.0": "1S0jrbULH9ZlXrDsxQPYms29REuYDNjDK",
    "32.0": "1MPWPvgmrbbzt-xmyra3mPcJFkXoYun3G",
}


@shrinked.command()
@click.option(
    "--version",
    "-v",
    "version",
    type=click.STRING,
    multiple=False,
    default="35.0",
    help="Version to install.",
)
def shrinked_download(version):
    """Download shrinked file, continue on break (it could take a few hours)."""
    if version not in versions:
        click.echo("Version " + version + " is not available as a skrinked DB.")
        return None
    filename = "./mysql/pfam" + version + ".sql.bz2"
    commands = [
        "mkdir -p mysql",
        "bash ./pfamserver/commands/db_shrinked_downloader.sh {id} {filename}",
    ]
    click.echo("\n".join(commands))
    # run(commands, id=versions[version], filename=filename)


@shrinked.command()
@click.option(
    "--version",
    "-v",
    "version",
    type=click.STRING,
    multiple=False,
    default="35.0",
    help="Version to install.",
)
def install(version):
    """Install shrinked file, continue on break (it could take a few hours)."""
    filename = "/home/pfamserver/stage/mysql/pfam" + version + ".sql.bz2"
    data_filename = "/home/pfamserver/stage/mysql/pfam" + version + ".sql"
    if not os.path.isfile(filename):
        click.echo("Version " + version + " wasn't downloaded as a shrinked DB.")
        return None
    commands = [
        "bzip2 -dk {filename}",
        'mysql -u root -h db --password="root" Pfam'
        + version.replace(".", "_")
        + " < {data_filename}",
    ]
    click.echo("\n".join(commands))
    # run(commands, filename=filename, data_filename=data_filename)
