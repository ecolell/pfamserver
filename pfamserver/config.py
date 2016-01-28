from application import app
from flask import send_from_directory
from database import db
# from autoupdate import update
from loader import init_db
import os


app.secret_key = 'some_secret_key!'
app.config['DEBUG'] = True
app.config['SECRET_KEY'] = app.secret_key
app.config.update(
    SQLALCHEMY_DATABASE_URI = 'mysql+pymysql://root:root@localhost:3306/pfamserver',
    #SQLALCHEMY_DATABASE_URI = 'sqlite:///test.db',
    #SQLALCHEMY_DATABASE_URI = 'postgresql://postgres:postgres@localhost:5432/pfamserver',
    SESSION_TYPE = 'sqlalchemy',
    SESSION_SQLALCHEMY = db
)

if not os.environ.get("WERKZEUG_RUN_MAIN"):
    init_db()
    # update()

@app.route('/favicon.ico')
def favicon():
    return send_from_directory(os.path.join(app.root_path, 'static/img'),
                               'fav.ico', mimetype='image/vnd.microsoft.icon')
