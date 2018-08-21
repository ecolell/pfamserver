from flask_restplus import Resource, abort

from pfamserver.api.v0 import api, schemas
from pfamserver.services import pfam_service
from pfamserver.extensions import cache
from flask import request

ns = api.namespace('pfams', decorators=[
    api.response(200, "success"),
    api.response(400, "not found")])


@ns.errorhandler(pfam_service.PfamServiceError)
def handle_root_exception(error):
    '''Return a custom message and 400 status code'''
    return {'message': error.message}, 400


def make_cache_key(*args, **kwargs):
    path = request.path
    args = str(request.args.items())
    return (path + args).encode('utf-8')


@ns.route('/<pfam>/sequence_descriptions')
class PfamAAPI(Resource):

    @ns.response(200, "response")
    @ns.doc('Obtain a sequence_description list from a pfam.')
    @cache.cached(timeout=3600, key_prefix=make_cache_key)
    @ns.expect(schemas.pfam_a_query)
    def get(self, pfam):
        kwargs = schemas.pfam_a_query.parse_args()
        with_pdb = kwargs['with_pdb']
        sequence_descriptions = pfam_service.get_sequence_descriptions_from_pfam(pfam, with_pdb)
        data = {'query': pfam,
                'with_pdb': with_pdb,
                'output': sequence_descriptions,
                'size': len(sequence_descriptions)}
        return data, 200


@ns.route('/<pfam>/stockholm')
class PfamAAPI(Resource):

    @ns.response(200, "response")
    @ns.doc('Obtain a sequence_description list from a pfam.')
    @cache.cached(timeout=3600, key_prefix=make_cache_key)
    @ns.expect(schemas.pfam_a_query)
    def get(self, pfam):
        kwargs = schemas.pfam_a_query.parse_args()
        with_pdb = kwargs['with_pdb']
        sequence_descriptions = pfam_service.get_sequence_descriptions_from_pfam(pfam, with_pdb)
        data = {'query': pfam,
                'with_pdb': with_pdb,
                'output': sequence_descriptions,
                'size': len(sequence_descriptions)}
        return data, 200
