from __future__ import unicode_literals
from pfamserver.database import db
from sqlalchemy import PrimaryKeyConstraint


class PfamARegFullSignificant(db.Model):
    __tablename__ = 'pfamA_reg_full_significant'

    pfamA_acc = db.Column(
        db.UnicodeText,
        index=True)
    pfamseq_acc = db.Column(
        db.UnicodeText,
        index=True)
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
