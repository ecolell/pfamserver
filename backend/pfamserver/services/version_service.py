import re

from flask import current_app


def version():
    database_uri = current_app.config.get("SQLALCHEMY_DATABASE_URI")
    ver = re.sub(r".*(Pfam\d+)_(\d+).*", r"\1.\2", database_uri)
    if "." not in ver:
        ver = f"Pfam{current_app.config.get('PFAM_VERSION')}"
    return ver


def hmmer_version():
    return current_app.config.get("HMMER_VERSION")
