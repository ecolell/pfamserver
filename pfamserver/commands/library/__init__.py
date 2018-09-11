# -*- coding: utf-8 -*-
from __future__ import absolute_import
from __future__ import print_function
from __future__ import unicode_literals

import click
click.disable_unicode_literals_warning = True
from pfamserver.commands.library.hmmer import hmmer as hmmer_command
from pfamserver.commands.library.pfamscan import pfamscan as pfamscan_command


@click.group()
def library():
    """Libraries commands"""
    pass


library.add_command(hmmer_command)
library.add_command(pfamscan_command)
