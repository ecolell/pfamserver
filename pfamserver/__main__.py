from application import app
import api
import admin
import web
import os


if __name__ == '__main__':
    host = os.getenv('HOST', '0.0.0.0')
    port = os.getenv('PORT', 5001)
    app.run(host=host, port=port, threaded=True)
