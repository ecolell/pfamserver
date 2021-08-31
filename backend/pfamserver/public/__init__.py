import os
from flask import Blueprint, send_from_directory, url_for
from flask import render_template

from pfamserver.extensions import cache

blueprint = Blueprint('', __name__, template_folder='templates')


@cache.memoize(timeout=86400)  # One day
@blueprint.route('/', defaults={'page': None})
def main_app(page):
    return render_template('index.html')


@blueprint.route('/favicon.ico')
def favicon():
    return send_from_directory(os.path.join(blueprint.root_path, 'static'),
                               'favicon.ico',
                               mimetype='image/png')
