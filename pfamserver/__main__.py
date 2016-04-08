import os

if __name__ == '__main__':
    from application import app
    from autoupdate import scheduler
    if os.environ.get("WERKZEUG_RUN_MAIN") == "true":
        scheduler.init_app(app)

        import api
        import admin
        import web
        from flask import request

        def shutdown_server():
            func = request.environ.get('werkzeug.server.shutdown')
            if func is None:
                raise RuntimeError('Not running with the Werkzeug Server')
            func()

    host = os.getenv('HOST', '0.0.0.0')
    port = int(os.getenv('PORT', '5001'))
    try:
        app.run(host=host, port=port, threaded=True)
    except KeyboardInterrupt:
        shutdown_server()
    if os.environ.get("WERKZEUG_RUN_MAIN") == "true":
        scheduler.finish_app(app)
