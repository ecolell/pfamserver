import sqlalchemy as sqla
from pfamserver.database import Model


class PfamARegFullSignificant(Model):
    __tablename__ = "pfamA_reg_full_significant"

    pfamA_acc = sqla.Column(sqla.String(7), index=True)
    pfamseq_acc = sqla.Column(sqla.String(10), index=True)
    seq_start = sqla.Column(sqla.Integer, index=True)
    seq_end = sqla.Column(sqla.Integer)
    in_full = sqla.Column(sqla.Integer, index=True)

    __table_args__ = (
        sqla.PrimaryKeyConstraint("pfamA_acc", "pfamseq_acc", "seq_start", "in_full"),
        {},
    )  # type: tuple
