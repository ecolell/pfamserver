from __future__ import unicode_literals
from pfamserver.database import db


class PfamAPfamseq(db.Model):
    pfamseq_id = db.Column(db.String(40)),
    pfamA_acc = db.Column(db.String(7)),
    pfamseq_acc = db.Column(db.String(10)),
    has_pdb = db.Column(db.Boolean, default=False)
