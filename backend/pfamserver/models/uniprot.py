import sqlalchemy as sqla
from pfamserver.database import Base


class Uniprot(Base):
    uniprot_acc = sqla.Column(sqla.String(10), primary_key=True)
    uniprot_id = sqla.Column(sqla.String(16), index=True, unique=True)
    description = sqla.Column(sqla.UnicodeText)
