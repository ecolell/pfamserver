from __future__ import unicode_literals

from flask_compress import Compress
from future.standard_library import install_aliases
install_aliases()

from flask_caching import Cache
from flask_sqlalchemy import SQLAlchemy
from flask_migrate import Migrate
from flask_bootstrap import Bootstrap
from flask_mail import Mail
from flask_webpack import Webpack
from flask_wtf.csrf import CsrfProtect
import requests
from six.moves.urllib.parse import quote
from base64 import b64decode
import json


bootstrap = Bootstrap()
cache = Cache(config={'CACHE_TYPE': 'simple'})
db = SQLAlchemy()
compress = Compress()
mail = Mail()
migrate = Migrate()
csrf = CsrfProtect()
webpack = Webpack()