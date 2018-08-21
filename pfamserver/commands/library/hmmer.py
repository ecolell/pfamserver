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
def install(version='3.2.1'):
    """Download hmmer sourcecode."""
    cmds = [
        'wget http://eddylab.org/software/hmmer/hmmer-{version}.tar.gz',
        'tar xvzf hmmer-{version}.tar.gz',
        'cd hmmer-{version}',
        './configure',
        'make'
    ]
    run(cmds, version=version)