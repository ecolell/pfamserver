from pfamserver.config.staging import env, StagingConfig


class ProductionConfig(StagingConfig):
    WEBPACK_ASSETS_URL = env.str('WEBPACK_ASSETS_URL', 'https://mistic2.leloir.org.ar:5001/static/dist/')

    @staticmethod
    def init_app(app):
        StagingConfig.init_app(app)