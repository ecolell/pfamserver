from typing import TYPE_CHECKING

import sqlalchemy as sqla
from pfamserver.database import Model
from sqlalchemy.dialects.mysql import INTEGER

if TYPE_CHECKING:
    from pfamserver.models.uniprot_reg_full import UniprotRegFull  # noqa: F401
    from pfamserver.models.pdb import Pdb  # noqa: F401


class PdbPfamAReg(Model):
    __tablename__ = "pdb_pfamA_reg"

    auto_uniprot_reg_full = sqla.Column(
        INTEGER(unsigned=True),
        sqla.ForeignKey("uniprot_reg_full.auto_uniprot_reg_full"),
    )
    uniprot_reg_full = sqla.orm.relationship(
        "UniprotRegFull", backref=sqla.orm.backref("pdbs")
    )
    pdb_id = sqla.Column(sqla.String(5), sqla.ForeignKey("pdb.pdb_id"))
    pdb = sqla.orm.relationship("Pdb", backref=sqla.orm.backref("pfams"))
    chain = sqla.Column(sqla.String(4), index=True)
    pfamA_acc = sqla.Column(sqla.String(7), index=True)
    pfamseq_acc = sqla.Column(sqla.String(10), index=True)
    pdb_res_start = sqla.Column(sqla.Integer, index=True)
    pdb_res_end = sqla.Column(sqla.Integer, index=True)

    __table_args__ = (
        sqla.PrimaryKeyConstraint("auto_uniprot_reg_full", "pdb_id", "chain"),
        {},
    )  # type: tuple
