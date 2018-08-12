from __future__ import unicode_literals
from pfamserver.database import db, Base


class PfamARegFullSignificant(Base):
    #auto_pfamA_reg_full = db.Column(
    #    db.Integer,
    #    primary_key=True)
    pfamA_acc = db.Column(
        db.UnicodeText,
        index=True,
        unique=True)
    pfamseq_acc = db.Column(
        db.UnicodeText,
        index=True,
        unique=True)
    seq_start = db.Column(
        db.Integer,
        index=True)
    seq_end = db.Column(
        db.Integer)
    in_full = db.Column(
        db.Integer,
        index=True)

    __table_args__ = (
        PrimaryKeyConstraint('pfamA_acc', 'pfamseq_acc', 'seq_start', 'in_full'),
        {},
    )