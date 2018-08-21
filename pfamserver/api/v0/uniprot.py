from flask_restplus import Resource, abort

from pfamserver.api.v0 import api, schemas
from pfamserver.services import uniprot_service
from pfamserver.extensions import cache


ns = api.namespace('uniprots', decorators=[
    api.response(200, "success"),
    api.response(400, "not found")])


@ns.errorhandler(uniprot_service.UniprotServiceError)
def handle_root_exception(error):
    '''Return a custom message and 400 status code'''
    return {'message': error.message}, 400


@ns.route('/<uniprot>/pfams')
class UniprotAPI(Resource):
    schema = schemas.UniprotSchema()

    @ns.response(200, "response")
    @ns.doc('Obtain a pfams list from a uniprot.')
    @cache.cached(timeout=3600)
    def get(self, uniprot):
        uniprot = uniprot_service.get_pfams_from_uniprot(uniprot)
        data, errors = self.schema.dump(uniprot)
        return data, 200
