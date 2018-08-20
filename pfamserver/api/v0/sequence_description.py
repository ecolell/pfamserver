from flask_restplus import Resource

from pfamserver.api.v0 import api, schemas
from pfamserver.services import pdb_service
from pfamserver.extensions import cache

ns = api.namespace('sequence_descriptions', decorators=[
    api.response(200, "success"),
    api.response(400, "not found")])


@ns.errorhandler(pdb_service.SequenceDescriptionServiceError)
def handle_root_exception(error):
    '''Return a custom message and 400 status code'''
    return {'message': error.message}, 400


@ns.route('/<uniprot_id>/<int:seq_start>-<int:seq_end>/pdbs')
class SequenceDescriptionAPI(Resource):
    schema = schemas.PdbPfamARegSchema()

    @ns.response(200, "response")
    @ns.doc('Obtain a pdb list from a sequence_description.')
    @cache.cached(timeout=3600)
    def get(self, uniprot_id, seq_start, seq_end):
        pdbs = pdb_service.get_pdbs_from_uniprot_pfam_a_reg(uniprot_id, seq_start, seq_end)
        output, errors = self.schema.dump(pdbs, many=True)
        data = {
            'query': {
                'uniprot_id': uniprot_id,
                'seq_start': seq_start,
                'seq_end': seq_end
            },
            'output': output
        }
        return data, 200