"""Pfamserver service."""
__version__ = "1.0.0"

import os
import sys

from flask import Flask
from pfamserver.config.development import DevelopmentConfig
from pfamserver.config.production import ProductionConfig
from pfamserver.config.staging import StagingConfig
from pfamserver.config.testing import TestingConfig
from pfamserver.extensions import cache, db

configs = {
    "development": DevelopmentConfig,
    "testing": TestingConfig,
    "staging": StagingConfig,
    "production": ProductionConfig,
    "default": DevelopmentConfig,
}


def create_app():
    """Build backend app."""
    config_name = os.getenv("FLASK_ENV", "default")
    config = configs.get(config_name)
    if config is None:
        print("Fail to initialize, unknown configuration {}".format(config_name))
        sys.exit(0)
    app = Flask(__name__)
    app.config.from_object(config)
    config.init_app(app)
    app.config["VERSION"] = __version__

    register_extensions(app)
    register_blueprints(app)

    register_middlewares(app)

    app.logger.setLevel(app.config["LOG_LEVEL"])
    app.logger.debug("Flask App created with '{}' config".format(config_name))

    return app


def register_blueprints(app):
    """Register blueprints into app."""
    from pfamserver import public
    from pfamserver.api.v0 import api_v0

    prefix = app.config.get("ROOT_URL_PREFIX", "")
    app.register_blueprint(public.blueprint)
    app.register_blueprint(api_v0, url_prefix=prefix + "/api/v0")

    return None


def register_extensions(app):
    """Register extensions into app."""
    db.init_app(app)
    app.db = db  # Should not be necessary, but was put here to support existing code.

    cache.init_app(app)


def register_middlewares(app):
    """Register middleware into app."""
    if app.config.get("PROXY_FIX"):
        from werkzeug.contrib.fixers import ProxyFix

        app.wsgi_app = ProxyFix(app.wsgi_app)

    if app.config.get("PROFILE", None):
        from werkzeug.contrib.profiler import ProfilerMiddleware

        app.wsgi_app = ProfilerMiddleware(app.wsgi_app, restrictions=[20])


def register_cli(app):
    """Register CLI commands into app."""
    from pfamserver.commands.library import library as library_command
    from pfamserver.commands.db import db as db_command

    app.cli.add_command(library_command)
    app.cli.add_command(db_command)
