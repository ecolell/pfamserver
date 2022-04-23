from pfamserver.database import db
from sqlalchemy import PrimaryKeyConstraint


class PfamAPfamseq(db.Model):
    __tablename__ = "pfamA_pfamseq"

    pfamseq_id = db.Column(db.String(40), index=True)
    pfamA_acc = db.Column(db.String(7), index=True)
    pfamseq_acc = db.Column(db.String(10), index=True)
    has_pdb = db.Column(db.Boolean, default=False)

    __table_args__ = (
        PrimaryKeyConstraint("pfamseq_id", "pfamA_acc", "pfamseq_acc"),
        {},
    )  # type: tuple
