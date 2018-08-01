from config.staging import env, StagingConfig


class ProductionConfig(StagingConfig):
    # Celery configuration
    CELERY_BROKER_URL = 'amqp://mistic:mistic@{:}:5672//'.format(StagingConfig.CELERY_SERVER_IP)
    CELERY_RESULT_BACKEND = 'amqp://mistic:mistic@{:}:5672//'.format(StagingConfig.CELERY_SERVER_IP)
    SQLALCHEMY_DATABASE_URI = env.str('SQLALCHEMY_DATABASE_URI', 'INVALID URL')

    @staticmethod
    def init_app(app):
        StagingConfig.init_app(app)