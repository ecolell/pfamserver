from flask_restx import Namespace, Resource
from pfamserver.api.v0 import schemas
from pfamserver.extensions import cache, make_cache_key
from pfamserver.services import sequence_service

ns = Namespace("protein_sequences")


@ns.errorhandler(sequence_service.SequenceServiceError)
def handle_root_exception(error):
    """Return a custom message and 400 status code"""
    return {"message": error.message}, 400


@ns.route("/<sequence>")
class ProteinSequenceAPI(Resource):
    schema = schemas.SequenceSchema()

    @ns.response(200, "response")
    @ns.doc("Obtain a pfams list from a uniprot.")
    @cache.cached(timeout=3600, make_cache_key=make_cache_key)
    def get(self, sequence):
        sequence = sequence.upper().strip()
        output = sequence_service.get_pfams_from_sequence(sequence)
        # mysql, errors = self.schema.dump(uniprot)
        data = {"query": sequence, "output": output}
        return data, 200
