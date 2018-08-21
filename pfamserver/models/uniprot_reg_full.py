from __future__ import unicode_literals
from pfamserver.database import db, Base
from sqlalchemy import PrimaryKeyConstraint


class UniprotRegFull(db.Model):
    auto_uniprot_reg_full = db.Column(
        db.Integer,
        primary_key=True)
    pfamA_acc = db.Column(
        db.UnicodeText(),
        db.ForeignKey('pfamA.pfamA_acc'))
    pfamA = db.relationship('PfamA', backref=db.backref('uniprots'))
    uniprot_acc = db.Column(
        db.UnicodeText(),
        db.ForeignKey('uniprot.uniprot_acc'))
    uniprot = db.relationship('Uniprot', backref=db.backref('pfams'))
    seq_start = db.Column(
        db.Integer,
        index=True)
    seq_end = db.Column(
        db.Integer)
    in_full = db.Column(
        db.Integer,
        index=True)

    __table_args__ = (
        PrimaryKeyConstraint('pfamA_acc', 'uniprot_acc', 'seq_start', 'in_full'),
        {},
    )
