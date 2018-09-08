from . import api
from marshmallow import Schema, fields
from flask_restplus import inputs


class UniprotRegFullSchema(Schema):
    pfamA_acc = fields.Str()
    description = fields.Str(attribute='pfamA.description')
    seq_start = fields.Int()
    seq_end = fields.Int()
    num_full = fields.Str(attribute='pfamA.num_full')


class UniprotSchema(Schema):
    uniprot_acc = fields.Str()
    uniprot_id = fields.Str()
    description = fields.Str()
    length = fields.Int()


class UniprotPfamsSchema(Schema):
    query = fields.Str(attribute='uniprot_id')
    output = fields.List(fields.Nested(UniprotRegFullSchema), attribute='pfams')


class PfamSchema(Schema):
    pfamA_acc = fields.Str()
    pfamA_id = fields.Str()
    description = fields.Str()
    num_full = fields.Int()


class PdbPfamARegSchema(Schema):
    pdb_id = fields.Str()
    chain = fields.Str()
    pdb_res_start = fields.Int()
    pdb_res_end = fields.Int()
    pfamA_acc = fields.Str()
    title = fields.Str(attribute='pdb.title')
    resolution = fields.Float(attribute='pdb.resolution')
    method = fields.Str(attribute='pdb.method')
    author = fields.Str(attribute='pdb.author')
    date = fields.Str(attribute='pdb.date')


pfam_a_query = api.parser()
pfam_a_query.add_argument('with_pdb', type=inputs.boolean, location='args', default=True)


class SequenceSchema(Schema):
    pass
