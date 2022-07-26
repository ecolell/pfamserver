import sqlalchemy as sqla
from pfamserver.database import Model


class PfamAPfamseq(Model):
    __tablename__ = "pfamA_pfamseq"

    pfamseq_id = sqla.Column(sqla.String(40), index=True)
    pfamA_acc = sqla.Column(sqla.String(7), index=True)
    pfamseq_acc = sqla.Column(sqla.String(10), index=True)
    has_pdb = sqla.Column(sqla.Boolean, default=False)

    __table_args__ = (
        sqla.PrimaryKeyConstraint("pfamseq_id", "pfamA_acc", "pfamseq_acc"),
        {},
    )  # type: tuple
