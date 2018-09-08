from flask_restplus import Resource

from pfamserver.api.v0 import api
from pfamserver.extensions import cache
from flask import request, current_app
import re


ns = api.namespace('version', decorators=[
    api.response(200, "success"),
    api.response(400, "not found")])


def make_cache_key(*args, **kwargs):
    path = request.path
    args = str(request.args.items())
    return (path + args).encode('utf-8')


@ns.route('')
class VersionAPI(Resource):

    @ns.response(200, "response")
    @ns.doc('Obtain the pfam database version.')
    @cache.cached(timeout=3600)
    def get(self):
        database_uri = current_app.config.get('SQLALCHEMY_DATABASE_URI')
        version = re.sub(r'.*(Pfam\d+)_(\d+).*', r'\1.\2', database_uri)
        data = {
            'version': version
        }
        return data, 200