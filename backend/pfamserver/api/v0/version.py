from flask_restplus import Resource

from pfamserver.api.v0 import api
from pfamserver.extensions import cache
from flask import request

from pfamserver.services import version_service

ns = api.namespace('version', decorators=[
    api.response(200, "success"),
    api.response(400, "not found")])


@ns.route('')
class VersionAPI(Resource):

    @ns.response(200, "response")
    @ns.doc('Obtain the pfam database version.')
    @cache.cached(timeout=3600)
    def get(self):
        data = {
            'version': version_service.version()
        }
        return data, 200
