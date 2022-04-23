from flask import Blueprint
from flask_restplus import Api
from pfamserver.api.v0.pfam import ns as pfam
from pfamserver.api.v0.protein_sequence import ns as protein_sequence
from pfamserver.api.v0.uniprot import ns as uniprot
from pfamserver.api.v0.version import ns as version
from pfamserver.extensions import csrf

api_v0 = Blueprint("api_v0", __name__)

api = Api(
    api_v0,
    version="0",
    title="pfamserver",
    description="API for the pfamserver project.",
    validate=True,
    decorators=[csrf.exempt],
)

api.add_namespace(pfam.ns, path="/pfams")
api.add_namespace(protein_sequence.ns, path="/protein_sequences")
api.add_namespace(uniprot.ns, path="/uniprots")
api.add_namespace(version.ns, path="/versions")
