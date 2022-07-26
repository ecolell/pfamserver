#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""PfamServer Wsgi booter."""
from __future__ import absolute_import
from __future__ import print_function
from __future__ import unicode_literals

from pfamserver import create_app

application = create_app()
app = application
