from flask_restplus import Namespace, Resource
from pfamserver.extensions import cache
from pfamserver.services import version_service

ns = Namespace("version")


@ns.route("")
class VersionAPI(Resource):
    @ns.response(200, "response")
    @ns.doc("Obtain the pfam database version.")
    @cache.cached(timeout=3600)
    def get(self):
        data = {"version": version_service.version()}
        return data, 200
