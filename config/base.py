import logging
import os
from datetime import timedelta
import logging

from environs import Env

env = Env()


class BaseConfig(object):
    # Serve bootstrap libs locally
    BOOTSTRAP_SERVE_LOCAL = True

    DEBUG = False
    LOG_SLOW_REQUESTS = False
    PROFILE = False
    PROFILE_SQL = False

    # LOG CONFIG
    LOG_LEVEL = logging.DEBUG
    LOG_FILE = '/tmp/flask_log.out'

    # Flask Toolbar
    DEBUG_TB_ENABLED = False
    DEBUG_TB_INTERCEPT_REDIRECTS = False

    # Flask Security
    SECRET_KEY = 'f49dc0a9de00d4e337e74b87bee8dbf806765ec73212ee7c9e6b1b5079a75eb7'
    SECURITY_PASSWORD_HASH = 'sha512_crypt'

    SECURITY_TRACKABLE = True
    SECURITY_RECOVERABLE = True
    SECURITY_SEND_REGISTER_EMAIL = False
    SECURITY_SEND_PASSWORD_CHANGE_EMAIL = False

    SQLALCHEMY_DATABASE_SCHEMA = 'public'
    SQLALCHEMY_ECHO = False  # True to log queries
    TABLE_CACHE_ENABLED = os.getenv('TABLE_CACHE_ENABLED', False)

    # Pagination setting
    ITEMS_PER_PAGE = 10

    # false not needed Flask-SQLAlchemy event system for now,
    # http://flask-sqlalchemy.pocoo.org/2.1/config/
    SQLALCHEMY_TRACK_MODIFICATIONS = False

    COLLECT_STATIC_ROOT = os.environ.get('COLLECT_STATIC_ROOT', 'public/static')

    BASE_PATH = os.path.dirname(os.path.realpath(__file__))
    UPLOAD_FOLDER = os.path.join(BASE_PATH, 'static/')

    WEBPACK_MANIFEST_PATH = './static/manifest.json'

    SEND_FILE_MAX_AGE_DEFAULT=3600

    SQLALCHEMY_DATABASE_URI = os.getenv('SQLALCHEMY_DATABASE_URI', 'mysql+pymysql://root:root@192.168.0.105:3306/Pfam32_0')
    PFAMSERVER_ROOT_PATH = os.getenv('PFAMSERVER_ROOT_PATH',
                                     'http://ftp.ebi.ac.uk/pub/databases/Pfam/releases/')

    # FLASK MAIL CONFIG
    MAIL_ENABLED = env.bool('MAIL_ENABLED', True)
    MAIL_SERVER = os.environ.get('MAIL_SERVER', 'smtp.gmail.com')
    MAIL_USE_TLS = env.bool('MAIL_USE_TLS', False)
    MAIL_USE_SSL = env.bool('MAIL_USE_SSL', True)
    MAIL_PORT = os.environ.get('MAIL_PORT', 465)
    MAIL_USERNAME = os.environ.get('MAIL_USERNAME', 'SET_THIS_ON_CONSOLE')
    MAIL_PASSWORD = os.environ.get('MAIL_PASSWORD', 'SET_THIS_ON_CONSOLE')
    MAIL_DEFAULT_SENDER = 'Pfamserver No Reply <{email}>'.format(email=MAIL_USERNAME)

    FEEDBACK_EMAIL = 'Feedback <mistic.leloir+pfamserver@gmail.com>'

    SENTRY_CONFIG = {}

    @staticmethod
    def init_app(app):
        pass