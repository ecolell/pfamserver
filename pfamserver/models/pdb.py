from __future__ import unicode_literals
from pfamserver.database import db


class Pdb(db.Model):
    pdb_id = db.Column(
        db.UnicodeText,
        primary_key=True)
    title = db.Column(db.UnicodeText)
    method = db.Column(db.UnicodeText)
    resolution = db.Column(db.Numeric)
    author = db.Column(db.UnicodeText)
    date = db.Column(db.UnicodeText)
