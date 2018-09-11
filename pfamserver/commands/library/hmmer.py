# -*- coding: utf-8 -*-
from __future__ import absolute_import
from __future__ import print_function
from __future__ import unicode_literals

import click
click.disable_unicode_literals_warning = True

import os


@click.group()
def hmmer():
    """HMMER commands"""
    pass


def run(cmds, **kwargs):
    cmd = ' && '.join(['({:})'.format(c) for c in cmds]).format(**kwargs)
    click.echo("Running:\n{}".format(cmd))
    os.system(cmd)


@hmmer.command()
@click.option('--version', '-v', 'version',
              type=click.STRING,
              multiple=False,
              default='3.2.1',
              help='Version to install.')
def install(version):
    """Download and compile hmmer sourcecode."""
    cmds = [
        'wget http://eddylab.org/software/hmmer/hmmer-{version}.tar.gz',
        'tar xvzf hmmer-{version}.tar.gz',
        'cd hmmer-{version} && ./configure',
        'cd hmmer-{version} && make',
        'ln -s hmmer-{version}/easel/miniapps/esl-afetch esl-afetch',
        'ln -s hmmer-{version}/src/hmmpres hmmpres',
        'ln -s hmmer-{version}/src/hmmscan hmmscan'
    ]
    run(cmds, version=version)


@hmmer.command()
@click.option('--version', '-v', 'version',
              type=click.STRING,
              multiple=False,
              default='31.0',
              help='Version to download.')
@click.option('--ftp',
              is_flag=True,
              help='Force to use ftp.')
def stockholm_index(version, ftp):
    """Download and index PfamA-full file."""
    protocol = 'ftp' if ftp else 'http'
    cmds = [
        'mkdir -p ./Pfam{version}',
        'wget -c {protocol}://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam{version}/Pfam-A.full.gz -O ./Pfam{version}/Pfam-A.full.gz',
        'gunzip -c ./Pfam{version}/Pfam-A.full.gz > ./Pfam{version}/Pfam-A.full',
        'sed -i -E "s/(#=GF AC   [A-Z0-9]+)\\.(.+)/\\1\\' + 'n#=GF DC   Revision: \\2/g" ./Pfam{version}/Pfam-A.full',
        './esl-afetch --index ./Pfam{version}/Pfam-A.full'
    ]
    run(cmds, version=version, protocol=protocol)
