from pfamserver.config.base import BaseConfig


class TestingConfig(BaseConfig):
    DEBUG = True
    TESTING = True

    @staticmethod
    def init_app(app):
        BaseConfig.init_app(app)
