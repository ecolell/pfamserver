from config.staging import env, StagingConfig


class ProductionConfig(StagingConfig):
    SQLALCHEMY_DATABASE_URI = env.str('SQLALCHEMY_DATABASE_URI', 'INVALID URL')

    @staticmethod
    def init_app(app):
        StagingConfig.init_app(app)