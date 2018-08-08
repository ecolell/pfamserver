from flask_restplus import fields
from . import api


pfam_a = api.model('Pfam', {
    'pfamA_acc': fields.String(attribute='UniprotRegFull.pfamA_acc'),
    'description': fields.String(attribute='PfamA.description'),
    'seq_start': fields.Integer(attribute='UniprotRegFull.seq_start'),
    'seq_end': fields.Integer(attribute='UniprotRegFull.seq_end'),
    'num_full': fields.String(attribute='PfamA.num_full')
})

pfams_from_uniprot = api.model('PfamsList', {
    'query': fields.String(attribute=lambda x: x[0].Uniprot.uniprot_id),
    'output': fields.List(fields.Nested(pfam_a), attribute=lambda x: x)
})