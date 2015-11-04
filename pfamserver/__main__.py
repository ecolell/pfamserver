from application import app
import config
import api
import web
from autoupdate import update


if __name__ == '__main__':
    update()
    app.run(host='0.0.0.0')
