from config.base import env, BaseConfig
import logging


class DevelopmentConfig(BaseConfig):
    DEBUG = True
    TESTING = False
    # BOOTSTRAP_SERVE_LOCAL = True
    MAIL_ENABLED = True

    @staticmethod
    def init_app(app):
        BaseConfig.init_app(app)

        # Handlers
        web_handler = logging.StreamHandler()
        web_fmt = logging.Formatter(
            "%(asctime)s %(levelname)s web %(name)s (%(process)d) %(message)s "
            "[in %(pathname)s:%(lineno)d]")
        web_handler.setFormatter(web_fmt)
        web_handler.setLevel(logging.DEBUG)

        worker_handler = logging.StreamHandler()
        worker_fmt = logging.Formatter(
            "%(asctime)s %(levelname)s worker %(name)s (%(process)d) %(message)s "
            "[in %(pathname)s:%(lineno)d]")
        worker_handler.setFormatter(worker_fmt)
        worker_handler.setLevel(logging.DEBUG)

        # Loggers
        app_logger = app.logger
        app_logger.setLevel(app.config['LOG_LEVEL'])
        app_logger.addHandler(web_handler)