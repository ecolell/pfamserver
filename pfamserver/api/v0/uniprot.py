from flask_restplus import Resource, abort

from pfamserver.api.v0 import api, schemas
from pfamserver.services import uniprot_service
from pfamserver.services import pdb_service
from pfamserver.extensions import cache


ns = api.namespace('uniprots', decorators=[
    api.response(200, "success"),
    api.response(400, "not found")])


@ns.errorhandler(uniprot_service.UniprotServiceError)
def handle_root_exception(error):
    '''Return a custom message and 400 status code'''
    return {'message': error.message}, 400


@ns.route('/<uniprot>')
class UniprotAPI(Resource):
    schema = schemas.UniprotSchema()

    @ns.response(200, "response")
    @ns.doc('Obtain the uniprot information.')
    @cache.cached(timeout=3600)
    def get(self, uniprot):
        uniprot = uniprot_service.get_uniprot(uniprot)
        data, errors = self.schema.dump(uniprot)
        return data, 200


@ns.route('/<uniprot>/pfams')
class UniprotPfamsAPI(Resource):
    schema = schemas.UniprotPfamsSchema()

    @ns.response(200, "response")
    @ns.doc('Obtain a pfams list from a uniprot.')
    @cache.cached(timeout=3600)
    def get(self, uniprot):
        uniprot = uniprot_service.get_pfams_from_uniprot(uniprot)
        data, errors = self.schema.dump(uniprot)
        return data, 200


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
