from pfamserver.config.base import BaseConfig


class StagingConfig(BaseConfig):
    @staticmethod
    def init_app(app):
        BaseConfig.init_app(app)
