from pfamserver.config.base import BaseConfig


class TestingConfig(BaseConfig):
    DEBUG = True
    TESTING = True
    # SQLALCHEMY_DATABASE_URI = env.str('SQLALCHEMY_DATABASE_URI', 'sqlite:///:memory:')

    @staticmethod
    def init_app(app):
        BaseConfig.init_app(app)
