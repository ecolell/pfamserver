from config.base import env, BaseConfig


class StagingConfig(BaseConfig):
    DEBUG = False
    TESTING = False

    @staticmethod
    def init_app(app):
        BaseConfig.init_app(app)