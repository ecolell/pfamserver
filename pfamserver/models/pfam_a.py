from __future__ import unicode_literals
from pfamserver.database import db, Base


class PfamA(Base):
    pfamA_acc = db.Column(
        db.UnicodeText,
        primary_key=True)
    pfamA_id = db.Column(
        db.UnicodeText,
        index=True,
        unique=True)
    description = db.Column(db.UnicodeText)
    num_full = db.Column(db.Integer)
