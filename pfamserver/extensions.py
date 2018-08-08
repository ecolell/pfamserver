from __future__ import unicode_literals

from future.standard_library import install_aliases

import logging
from flask_caching import Cache
from flask_collect import Collect
from flask_compress import Compress
from flask_wtf.csrf import CSRFProtect
from flask_sqlalchemy import SQLAlchemy
from flask_security import Security
from flask_webpack import Webpack
from raven.contrib.flask import Sentry

install_aliases()

cache = Cache(config={'CACHE_TYPE': 'simple'})
collect = Collect()
compress = Compress()
csrf = CSRFProtect()
db = SQLAlchemy()
sentry = Sentry(logging=True, level=logging.ERROR)
webpack = Webpack()