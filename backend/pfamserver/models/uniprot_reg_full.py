from typing import TYPE_CHECKING

import sqlalchemy as sqla
from pfamserver.database import Model
from sqlalchemy.dialects.mysql import INTEGER

if TYPE_CHECKING:
    from pfamserver.models.pfam_a import PfamA  # noqa: F401
    from pfamserver.models.uniprot import Uniprot  # noqa: F401


class UniprotRegFull(Model):
    auto_uniprot_reg_full = sqla.Column(INTEGER(unsigned=True), primary_key=True)
    pfamA_acc = sqla.Column(sqla.String(7), sqla.ForeignKey("pfamA.pfamA_acc"))
    pfamA = sqla.orm.relationship("PfamA", backref=sqla.orm.backref("uniprots"))
    uniprot_acc = sqla.Column(sqla.String(10), sqla.ForeignKey("uniprot.uniprot_acc"))
    uniprot = sqla.orm.relationship("Uniprot", backref=sqla.orm.backref("pfams"))
    seq_start = sqla.Column(sqla.Integer, index=True)
    seq_end = sqla.Column(sqla.Integer)
    in_full = sqla.Column(sqla.Integer, index=True)

    __table_args__ = (
        sqla.PrimaryKeyConstraint("pfamA_acc", "uniprot_acc", "seq_start", "in_full"),
        {},
    )  # type: tuple
