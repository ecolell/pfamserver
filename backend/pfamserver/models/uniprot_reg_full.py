from pfamserver.database import db
from sqlalchemy import PrimaryKeyConstraint
from sqlalchemy.dialects.mysql import INTEGER


class UniprotRegFull(db.Model):
    auto_uniprot_reg_full = db.Column(INTEGER(unsigned=True), primary_key=True)
    pfamA_acc = db.Column(db.String(7), db.ForeignKey("pfamA.pfamA_acc"))
    pfamA = db.relationship("PfamA", backref=db.backref("uniprots"))
    uniprot_acc = db.Column(db.String(10), db.ForeignKey("uniprot.uniprot_acc"))
    uniprot = db.relationship("Uniprot", backref=db.backref("pfams"))
    seq_start = db.Column(db.Integer, index=True)
    seq_end = db.Column(db.Integer)
    in_full = db.Column(db.Integer, index=True)

    __table_args__ = (
        PrimaryKeyConstraint("pfamA_acc", "uniprot_acc", "seq_start", "in_full"),
        {},
    )  # type: tuple
