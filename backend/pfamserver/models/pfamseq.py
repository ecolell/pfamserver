import sqlalchemy as sqla
from pfamserver.database import Base


class Pfamseq(Base):
    pfamseq_acc = sqla.Column(sqla.String(10), primary_key=True)
    pfamseq_id = sqla.Column(sqla.String(16), index=True)
