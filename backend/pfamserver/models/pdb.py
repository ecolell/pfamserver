import sqlalchemy as sqla
from pfamserver.database import Model


class Pdb(Model):
    pdb_id = sqla.Column(sqla.String(5), primary_key=True)
    title = sqla.Column(sqla.UnicodeText)
    method = sqla.Column(sqla.UnicodeText)
    resolution = sqla.Column(sqla.Numeric)
    author = sqla.Column(sqla.UnicodeText)
    date = sqla.Column(sqla.UnicodeText)
