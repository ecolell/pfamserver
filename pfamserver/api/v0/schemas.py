from marshmallow import Schema, fields


class PfamASchema(Schema):
    pfamA_acc = fields.Str()
    description = fields.Str(attribute='pfamA.description')
    seq_start = fields.Int()
    seq_end =  fields.Int()
    num_full = fields.Str(attribute='pfamA.num_full')