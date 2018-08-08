from config.base import env, BaseConfig


class StagingConfig(BaseConfig):
    DEBUG = False
    TESTING = False
    SQLALCHEMY_DATABASE_URI = env.str('SQLALCHEMY_DATABASE_URI',
                                      'postgres:///pfamserver:password@localhost:5432/pfamserverstage')

    @staticmethod
    def init_app(app):
        BaseConfig.init_app(app)