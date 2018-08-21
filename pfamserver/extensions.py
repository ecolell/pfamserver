from __future__ import unicode_literals

import logging
from flask_caching import Cache
from flask_collect import Collect
from flask_compress import Compress
from flask_wtf.csrf import CSRFProtect
from flask_sqlalchemy import SQLAlchemy
from flask_webpack import Webpack
from raven.contrib.flask import Sentry

cache = Cache(config={'CACHE_TYPE': 'simple'})
collect = Collect()
compress = Compress()
csrf = CSRFProtect()
db = SQLAlchemy()  # session_options={"autoflush": False})
sentry = Sentry(logging=True, level=logging.ERROR)
webpack = Webpack()
