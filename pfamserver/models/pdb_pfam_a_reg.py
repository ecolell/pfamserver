from __future__ import unicode_literals
from pfamserver.database import db
from sqlalchemy import PrimaryKeyConstraint


class PdbPfamAReg(db.Model):
    __tablename__ = 'pdb_pfamA_reg'

    # auto_pfamA_reg_full = db.Column(
    #    db.Integer,
    #    primary_key=True)
    auto_uniprot_reg_full = db.Column(
        db.UnicodeText,
        index=True)
    pdb_id = db.Column(
        db.UnicodeText,
        index=True)
    chain = db.Column(
        db.UnicodeText,
        index=True)
    pfamA_acc = db.Column(
        db.UnicodeText,
        index=True)
    pfamseq_acc = db.Column(
        db.UnicodeText,
        index=True
    )
    pdb_res_start = db.Column(
        db.Integer,
        index=True)
    pdb_res_end = db.Column(
        db.Integer,
        index=True)

    __table_args__ = (
        PrimaryKeyConstraint('auto_uniprot_reg_full', 'pdb_id'),
        {},
    )