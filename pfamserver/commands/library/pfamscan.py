# -*- coding: utf-8 -*-
from __future__ import absolute_import
from __future__ import print_function
from __future__ import unicode_literals

import click
from flask.cli import with_appcontext
from flask import current_app


click.disable_unicode_literals_warning = True

import os


@click.group()
def pfamscan():
    """PfamScan commands"""
    pass


def run(cmds, **kwargs):
    cmd = ' && '.join(['({:})'.format(c) for c in cmds]).format(**kwargs)
    click.echo("Running:\n{}".format(cmd))
    os.system(cmd)


@pfamscan.command()
def install():
    """Download and compile pfamscan sourcecode."""
    cmds = [
        'wget -c http://ftp.ebi.ac.uk/pub/databases/Pfam/Tools/PfamScan.tar.gz',
        'tar xvzf PfamScan.tar.gz',
        'sudo perl - MCPAN - e"install Moose"'
    ]
    run(cmds)


@pfamscan.command()
@click.option('--version', '-v', 'version',
              type=click.STRING,
              multiple=False,
              default='31.0',
              help='Version to download.')
@click.option('--ftp',
              is_flag=True,
              help='Force to use ftp.')
def index(version, ftp):
    """Download PfamA-full file."""
    protocol = 'ftp' if ftp else 'http'
    cmds = [
        'wget -c {protocol}://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam{version}/Pfam-A.hmm.gz -O ./Pfam{version}/Pfam-A.hmm.gz',
        'wget -c {protocol}://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam{version}/Pfam-A.hmm.dat.gz -O ./Pfam{version}/Pfam-A.hmm.dat.gz',
        'gunzip -c ./Pfam{version}/Pfam-A.hmm.gz > ./Pfam{version}/Pfam-A.hmm',
        'gunzip -c ./Pfam{version}/Pfam-A.hmm.dat.gz > ./Pfam{version}/Pfam-A.hmm.dat',
        './hmmpress ./Pfam{version}/Pfam-A.hmm',
        'mkdir -p tmp',

    ]
    run(cmds, version=version, protocol=protocol)
