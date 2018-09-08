from flask_restplus import Resource, abort

from pfamserver.api.v0 import api, schemas
from pfamserver.services import sequence_service
from pfamserver.extensions import cache


ns = api.namespace('protein_sequences', decorators=[
    api.response(200, "success"),
    api.response(400, "not found")])


@ns.errorhandler(sequence_service.SequenceServiceError)
def handle_root_exception(error):
    '''Return a custom message and 400 status code'''
    return {'message': error.message}, 400


@ns.route('/<sequence>')
class UniprotAPI(Resource):
    schema = schemas.SequenceSchema()

    @ns.response(200, "response")
    @ns.doc('Obtain a pfams list from a uniprot.')
    @cache.cached(timeout=3600)
    def get(self, sequence):
        uniprot = sequence_service.get_pfams_from_sequence(sequence)
        data, errors = self.schema.dump(uniprot)
        return data, 200
