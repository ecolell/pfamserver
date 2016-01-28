from application import app
import config
import admin
import api
import web
import os


if __name__ == '__main__':
    host = os.getenv('IP', '0.0.0.0')
    port = os.getenv('PORT', 5001)
    app.run(host=host, port=port, threaded=True)
