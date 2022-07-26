import sqlalchemy as sqla
from pfamserver.database import Base


class PfamA(Base):
    pfamA_acc = sqla.Column(sqla.String(7), primary_key=True)
    pfamA_id = sqla.Column(sqla.String(16), index=True, unique=True)
    description = sqla.Column(sqla.UnicodeText)
    num_full = sqla.Column(sqla.Integer)
