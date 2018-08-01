from config.base import env, BaseConfig


class StagingConfig(BaseConfig):
    DEBUG = False
    TESTING = False
    # Celery configuration
    CELERY_BROKER_URL = 'amqp://mistic:mistic@{:}:5672//'.format(BaseConfig.CELERY_SERVER_IP)
    CELERY_RESULT_BACKEND = 'amqp://mistic:mistic@{:}:5672//'.format(BaseConfig.CELERY_SERVER_IP)
    SQLALCHEMY_DATABASE_URI = env.str('SQLALCHEMY_DATABASE_URI',
                                      'postgres:///pfamserver:password@localhost:5432/pfamserverstage')

    @staticmethod
    def init_app(app):
        BaseConfig.init_app(app)