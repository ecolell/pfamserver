import os
from multiprocessing import Process


def create_server():
    from application import app
    from autoupdate import scheduler
    debugging = app.config['DEBUG']
    run_once = os.environ.get("WERKZEUG_RUN_MAIN") == "true" or not debugging


    def shutdown_server():
        func = request.environ.get('werkzeug.server.shutdown')
        if func is None:
            raise RuntimeError('Not running with the Werkzeug Server')
        func()
    app.shutdown = shutdown_server

    if run_once:
        scheduler.init_app(app)

        import api
        import admin
        import web
        from flask import request

    host = os.getenv('HOST', '0.0.0.0')
    port = int(os.getenv('PORT', '5001'))
    try:
        app.run(host=host, port=port, threaded=True)
    except KeyboardInterrupt, e:
        if debugging:
            shutdown_server()
        raise KeyboardInterrupt(e)
    if run_once:
        scheduler.finish_app(app)


if __name__ == '__main__':
    create_server()
