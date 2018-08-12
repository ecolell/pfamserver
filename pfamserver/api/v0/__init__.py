from flask import Blueprint
from flask_restplus import Api
from pfamserver.extensions import csrf

api_v0 = Blueprint('api_v0', __name__)

api = Api(
    api_v0,
    version='0',
    title='pfamserver',
    description='API for the pfamserver project.',
    validate=True,
    decorators=[csrf.exempt]
)


from pfamserver.api.v0.uniprot import ns as uniprot
#from pfamserver.api.v0.msa_settings import ns as msa_settings
#from pfamserver.api.v0.pdb_settings import ns as pdb_settings
#from pfamserver.api.v0.algorithms import ns as algorithms
#from pfamserver.api.v0.algorithm_settings import ns as algorithm_settings
#from pfamserver.api.v0.notifications import ns as notifications
