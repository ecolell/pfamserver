from base64 import b64encode
from zlib import compress

from flask import request
from flask_restx import Namespace, Resource
from pfamserver.api.v0 import schemas
from pfamserver.extensions import cache
from pfamserver.services import pfam_service
from webargs.flaskparser import use_kwargs

ns = Namespace("pfams")


@ns.errorhandler(pfam_service.PfamServiceError)
def handle_root_exception(error):
    """Return a custom message and 400 status code"""
    return {"message": error.message}, 400


def make_cache_key(*args, **kwargs):
    path = request.path
    args = str(request.args.items())
    return (path + args).encode("utf-8")


@ns.route("/<pfam>")
class PfamAAPI(Resource):
    schema = schemas.PfamSchema()

    @ns.response(200, "response")
    @ns.doc("Obtain the pfam information.")
    @cache.cached(timeout=3600)
    def get(self, pfam):
        pfam = pfam_service.get_pfam(pfam)
        data = self.schema.dump(pfam)
        return data, 200


@ns.route("/<pfam>/sequence_descriptions")
class PfamASequenceDescriptionsAPI(Resource):
    @ns.response(200, "response")
    @ns.doc("Obtain a sequence_description list from a pfam.")
    @cache.cached(timeout=3600, key_prefix=make_cache_key)
    @use_kwargs(schemas.PfamAQuery, location="args")
    def get(self, pfam):
        kwargs = schemas.pfam_a_query.parse_args()
        with_pdb = kwargs["with_pdb"]
        sequence_descriptions = pfam_service.get_sequence_descriptions_from_pfam(
            pfam, with_pdb
        )
        data = {
            "query": pfam,
            "with_pdb": with_pdb,
            "output": sequence_descriptions,
            "size": len(sequence_descriptions),
        }
        return data, 200


@ns.route("/<pfam>/stockholm")
class PfamAStockholmAPI(Resource):
    @ns.response(200, "response")
    @ns.doc("Obtain a sequence_description list from a pfam.")
    @cache.cached(timeout=3600, key_prefix=make_cache_key)
    def get(self, pfam):
        stockholm = pfam_service.get_stockholm_from_pfam(pfam)
        data = {"query": pfam, "output": b64encode(compress(stockholm))}
        return data, 200