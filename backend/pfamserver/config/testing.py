import os

from pfamserver.config.base import BaseConfig


class TestingConfig(BaseConfig):
    DEBUG = True
    TESTING = True

    CACHE_TYPE = os.environ.get("CACHE_TYPE", "NullCache")

    @staticmethod
    def init_app(app):
        BaseConfig.init_app(app)
