from pfamserver.config.base import env, BaseConfig


class StagingConfig(BaseConfig):
    WEBPACK_ASSETS_URL = env.str('WEBPACK_ASSETS_URL', 'https://mistic2.herokuapp.com:5001/static/dist/')

    @staticmethod
    def init_app(app):
        BaseConfig.init_app(app)