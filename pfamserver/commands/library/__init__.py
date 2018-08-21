# -*- coding: utf-8 -*-
from __future__ import absolute_import
from __future__ import print_function
from __future__ import unicode_literals

import click
click.disable_unicode_literals_warning = True
from pfamserver.commands.library.hmmer import hmmer as hmmer_command


@click.group()
def library():
    """Data commands"""
    pass


library.add_command(hmmer_command)
