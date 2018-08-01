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
    SECRET_KEY = 'a6a5c4dce409af370a1c94ab07c447a365dadfe99a8d5df635a556f9'
    OAUTH2_PROVIDER_TOKEN_EXPIRES_IN = 3600 * 24 * 365 * 10  # 10 years in seconds
    SECURITY_PASSWORD_HASH = 'sha512_crypt'

    SECURITY_TRACKABLE = True
    SECURITY_RECOVERABLE = True
    SECURITY_SEND_REGISTER_EMAIL = False
    SECURITY_SEND_PASSWORD_CHANGE_EMAIL = False

    SQLALCHEMY_DATABASE_SCHEMA = 'public'
    SQLALCHEMY_ECHO = False  # True to log queries

    # Pagination setting
    ITEMS_PER_PAGE = 10

    # false not needed Flask-SQLAlchemy event system for now,
    # http://flask-sqlalchemy.pocoo.org/2.1/config/
    SQLALCHEMY_TRACK_MODIFICATIONS = False

    COLLECT_STATIC_ROOT = os.environ.get('COLLECT_STATIC_ROOT', 'public/static')

    BASE_PATH = os.path.dirname(os.path.realpath(__file__))
    UPLOAD_FOLDER = os.path.join(BASE_PATH, 'static/')

    WEBPACK_MANIFEST_PATH = './static/manifest.json'

    CELERYBEAT_SCHEDULE = {
        'detect-new-version-every-1-hour': {
            'task': 'register_new_versions',
            'schedule': timedelta(hours=1)
        },
    #    'notifier-jobrequest-every-10-seconds': {
    #        'task': 'job_request.notifier',
    #        'schedule': timedelta(seconds=10)
    #    },
    #    'clean-jobrequest-every-2-hours': {
    #        'task': 'job_request.clean_open_jobs',
    #        'schedule': timedelta(hours=2)
    #    },
    #    'clean-jobrequest-every-2-hours': {
    #        'task': 'job_request.clean_finished_jobs',
    #        'schedule': timedelta(days=30)
    #    },
    }

    CELERY_CORES = 3

    SEND_FILE_MAX_AGE_DEFAULT=3600

    SQLALCHEMY_DATABASE_URI = os.getenv('SQLALCHEMY_DATABASE_URI')
    PFAMSERVER_ROOT_PATH = os.getenv('PFAMSERVER_ROOT_PATH',
                                     'http://ftp.ebi.ac.uk/pub/databases/Pfam/releases/')


    # Celery configuration
    CELERY_SERVER_IP = env.str('CELERY_SERVER_IP', 'localhost')
    CELERY_BROKER_URL = 'amqp://guest:guest@{:}:5672//'.format(CELERY_SERVER_IP)

    @staticmethod
    def init_app(app):
        pass