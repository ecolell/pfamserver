from flask import current_app
import re


def version():
    database_uri = current_app.config.get('SQLALCHEMY_DATABASE_URI')
    return re.sub(r'.*(Pfam\d+)_(\d+).*', r'\1.\2', database_uri)
