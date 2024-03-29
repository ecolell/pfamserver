import logging
import os


class BaseConfig:
    # Serve bootstrap libs locally
    BOOTSTRAP_SERVE_LOCAL = True

    DEBUG = False
    LOG_SLOW_REQUESTS = False
    PROFILE = False
    PROFILE_SQL = False

    # LOG CONFIG
    LOG_LEVEL = logging.DEBUG

    # Flask Toolbar
    DEBUG_TB_ENABLED = False
    DEBUG_TB_INTERCEPT_REDIRECTS = False

    SQLALCHEMY_DATABASE_SCHEMA = "public"
    SQLALCHEMY_ECHO = False  # True to log queries
    TABLE_CACHE_ENABLED = os.getenv("TABLE_CACHE_ENABLED", True)

    # Caching
    CACHE_TYPE = os.environ.get("CACHE_TYPE", "RedisCache")
    CACHE_REDIS_HOST = os.environ.get("CACHE_REDIS_HOST", "redis-cache")
    CACHE_REDIS_PORT = os.environ.get("CACHE_REDIS_PORT", 6379)
    CACHE_REDIS_DB = os.environ.get("CACHE_REDIS_DB", 0)
    CACHE_REDIS_URL = os.environ.get("CACHE_REDIS_URL", "redis://redis-cache:6379/0")
    CACHE_DEFAULT_TIMEOUT = os.environ.get("CACHE_DEFAULT_TIMEOUT", 500)

    # Pagination setting
    ITEMS_PER_PAGE = 10

    # false not needed Flask-SQLAlchemy event system for now,
    # http://flask-sqlalchemy.pocoo.org/2.1/config/
    SQLALCHEMY_TRACK_MODIFICATIONS = False

    COLLECT_STATIC_ROOT = os.environ.get("COLLECT_STATIC_ROOT", "public/static")

    BASE_PATH = os.path.dirname(os.path.realpath(__file__))
    UPLOAD_FOLDER = os.path.join(BASE_PATH, "static/")

    SEND_FILE_MAX_AGE_DEFAULT = 3600

    PFAM_VERSION = os.getenv("PFAM_VERSION", "35.0")
    HMMER_VERSION = os.getenv("HMMER_VERSION", "3.2.1")
    SQLALCHEMY_DATABASE_URI = os.getenv(
        "SQLALCHEMY_DATABASE_URI", "mysql+pymysql://root:root@db:3306/Pfam35_0"
    )
    PFAMSERVER_ROOT_PATH = os.getenv(
        "PFAMSERVER_ROOT_PATH", "http://ftp.ebi.ac.uk/pub/databases/Pfam/releases/"
    )

    @staticmethod
    def init_app(app):
        pass
