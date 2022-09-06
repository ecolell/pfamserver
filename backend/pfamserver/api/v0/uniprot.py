from flask_restx import Namespace, Resource
from pfamserver.api.v0 import schemas
from pfamserver.extensions import cache, make_cache_key
from pfamserver.services import pdb_service, uniprot_service

ns = Namespace("uniprots")


@ns.errorhandler(uniprot_service.UniprotServiceError)
def handle_root_exception(error):
    """Return a custom message and 400 status code"""
    return {"message": error.message}, 400


@ns.route("/<uniprot>")
class UniprotAPI(Resource):
    schema = schemas.UniprotSchema()

    @ns.response(200, "response")
    @ns.doc("Obtain the uniprot information.")
    @cache.cached(timeout=3600, key_prefix=make_cache_key)
    def get(self, uniprot):
        uniprot = uniprot_service.get_uniprot(uniprot)
        data = self.schema.dump(uniprot)
        return data, 200


@ns.route("/<uniprot>/pfams")
class UniprotPfamsAPI(Resource):
    schema = schemas.UniprotPfamsSchema()

    @ns.response(200, "response")
    @ns.doc("Obtain a pfams list from a uniprot.")
    @cache.cached(timeout=3600, key_prefix=make_cache_key)
    def get(self, uniprot):
        uniprot = uniprot_service.get_pfams_from_uniprot(uniprot)
        data = self.schema.dump(uniprot)
        return data, 200


@ns.route("/<uniprot_id>/<int:seq_start>-<int:seq_end>/pdbs")
class SequenceDescriptionAPI(Resource):
    schema = schemas.PdbPfamARegSchema()

    @ns.response(200, "response")
    @ns.doc("Obtain a pdb list from a sequence_description.")
    @cache.cached(timeout=3600, key_prefix=make_cache_key)
    def get(self, uniprot_id, seq_start, seq_end):
        pdbs = pdb_service.get_pdbs_from_uniprot_pfam_a_reg(
            uniprot_id, seq_start, seq_end
        )
        output = self.schema.dump(pdbs, many=True)
        data = {
            "query": {
                "uniprot_id": uniprot_id,
                "seq_start": seq_start,
                "seq_end": seq_end,
            },
            "output": output,
        }
        return data, 200
