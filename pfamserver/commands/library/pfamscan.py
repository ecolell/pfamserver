# -*- coding: utf-8 -*-
from __future__ import absolute_import
from __future__ import print_function
from __future__ import unicode_literals

import click
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
@click.option('--version', '-v', 'version',
              type=click.STRING,
              multiple=False,
              default='3.2.1',
              help='Version to install.')
def install(version):
    """Download and compile pfamscan sourcecode."""
    cmds = [
        'wget http://eddylab.org/software/pfamscan/pfamscan-{version}.tar.gz',
        'tar xvzf pfamscan-{version}.tar.gz',
        'cd pfamscan-{version} && ./configure',
        'cd pfamscan-{version} && make',
        'ln -s pfamscan-{version}/easel/miniapps/esl-afetch esl-afetch'
    ]
    run(cmds, version=version)


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
        'wget -c {protocol}://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam{version}/Pfam-A.hmm.gz',
        'wget -c {protocol}://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam{version}/Pfam-A.hmm.dat.gz',
        'gunzip Pfam-A.hmm*.gz',
    ]
    run(cmds, version=version, protocol=protocol)
