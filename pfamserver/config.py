from application import app
from flask import send_from_directory
import os


app.secret_key = 'some_secret_key!'
app.config['DEBUG'] = True


@app.route('/favicon.ico')
def favicon():
    return send_from_directory(os.path.join(app.root_path, 'static/img'),
                               'fav.ico', mimetype='image/vnd.microsoft.icon')