import os

from flask import Blueprint, render_template, send_from_directory
from pfamserver.extensions import cache

blueprint = Blueprint("", __name__, template_folder="templates")


@cache.memoize(timeout=86400)  # One day
@blueprint.route("/")
def main_app():
    return render_template("index.html")


@blueprint.route("/favicon.ico")
def favicon():
    return send_from_directory(
        os.path.join(blueprint.root_path, "static"), "favicon.ico", mimetype="image/png"
    )
