from pfamserver.config.staging import StagingConfig


class ProductionConfig(StagingConfig):
    @staticmethod
    def init_app(app):
        StagingConfig.init_app(app)
