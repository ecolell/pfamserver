from pfamserver.database import db
from sqlalchemy import PrimaryKeyConstraint
from sqlalchemy.dialects.mysql import INTEGER


class PdbPfamAReg(db.Model):
    __tablename__ = "pdb_pfamA_reg"

    # auto_pfamA_reg_full = db.Column(
    #    db.Integer,
    #    primary_key=True)
    auto_uniprot_reg_full = db.Column(
        INTEGER(unsigned=True), db.ForeignKey("uniprot_reg_full.auto_uniprot_reg_full")
    )
    uniprot_reg_full = db.relationship("UniprotRegFull", backref=db.backref("pdbs"))
    pdb_id = db.Column(db.String(5), db.ForeignKey("pdb.pdb_id"))
    pdb = db.relationship("Pdb", backref=db.backref("pfams"))
    chain = db.Column(db.String(4), index=True)
    pfamA_acc = db.Column(db.String(7), index=True)
    pfamseq_acc = db.Column(db.String(10), index=True)
    pdb_res_start = db.Column(db.Integer, index=True)
    pdb_res_end = db.Column(db.Integer, index=True)

    __table_args__ = (
        PrimaryKeyConstraint("auto_uniprot_reg_full", "pdb_id"),
        {},
    )  # type: tuple
