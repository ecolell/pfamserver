from . import api
from marshmallow import Schema, fields
from flask_restplus import reqparse
from flask_restplus import inputs


class UniprotRegFullSchema(Schema):
    pfamA_acc = fields.Str()
    description = fields.Str(attribute='pfamA.description')
    seq_start = fields.Int()
    seq_end = fields.Int()
    num_full = fields.Str(attribute='pfamA.num_full')


class UniprotSchema(Schema):
    query = fields.Str(attribute='uniprot_id')
    output = fields.List(fields.Nested(UniprotRegFullSchema), attribute='pfams')


pfam_a_query = reqparse.RequestParser()
pfam_a_query.add_argument('with_pdb', type=inputs.boolean, location='args', default='true')