import os

if __name__ == '__main__':
    from application import app
    if os.environ.get("WERKZEUG_RUN_MAIN") == "true":
        from autoupdate import scheduler
        scheduler.init_app(app)

        import api
        import admin
        import web
    host = os.getenv('HOST', '0.0.0.0')
    port = os.getenv('PORT', 5001)
    app.run(host=host, port=port, threaded=True)
