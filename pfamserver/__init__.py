__version__ = '1.0.0'

import sys
import os
from pfamserver import static
from pfamserver.extensions import db, cache, collect, compress, webpack, sentry
from flask import Flask
from flask_restplus import apidoc

from config.development import DevelopmentConfig
from config.staging import StagingConfig
from config.production import ProductionConfig
from config.testing import TestingConfig
from pfamserver.exceptions import SentryIgnoredError

configs = {
    'development': DevelopmentConfig,
    'testing': TestingConfig,
    'staging': StagingConfig,
    'production': ProductionConfig,
    'default': DevelopmentConfig
}


def create_app():
    config_name = os.getenv('FLASK_ENV', 'default')
    config = configs.get(config_name)
    if config is None:
        print('Fail to initialize, unknown configuration {}'.format(config_name))
        sys.exit(0)
    app = Flask(__name__)
    app.config.from_object(config)
    config.init_app(app)
    app.config['VERSION'] = __version__

    app.config['SENTRY_CONFIG']['release'] = __version__
    app.config['SENTRY_CONFIG']['ignore_exceptions'] = [
        SentryIgnoredError
    ]

    register_extensions(app)
    register_blueprints(app)
    register_cli(app)

    register_middlewares(app)

    app.logger.setLevel(app.config['LOG_LEVEL'])
    app.logger.debug("Flask App created with '{0}' config".format(config_name))

    return app


def register_blueprints(app):
    "Registers blueprints into app"
    from pfamserver.api.v0 import api_v0

    prefix = app.config.get('ROOT_URL_PREFIX', '')
    app.register_blueprint(static.blueprint)
    app.register_blueprint(api_v0, url_prefix=prefix + '/api/v0')

    if app.config.get('DEBUG'):
        app.register_blueprint(apidoc.apidoc, url_prefix=prefix)

    return None


def register_extensions(app):
    db.init_app(app)
    app.db = db  # Should not be necessary, but was put here to support existing code.

    cache.init_app(app)
    collect.init_app(app)
    compress.init_app(app)
    sentry.init_app(app)
    webpack.init_app(app)


def register_middlewares(app):
    """Middleware registration"""

    if app.config.get('PROXY_FIX'):
        from werkzeug.contrib.fixers import ProxyFix
        app.wsgi_app = ProxyFix(app.wsgi_app)

    if app.config.get('PROFILE', None):
        from werkzeug.contrib.profiler import ProfilerMiddleware
        app.wsgi_app = ProfilerMiddleware(app.wsgi_app, restrictions=[20])

    if app.config.get('PROFILE_SQL', None):
        from sqltap.wsgi import SQLTapMiddleware
        app.wsgi_app = SQLTapMiddleware(app.wsgi_app)


def register_cli(app):
    from pfamserver.commands.library import library as library_command
    from pfamserver.commands.data import data as data_command

    app.cli.add_command(library_command)
    app.cli.add_command(data_command)
