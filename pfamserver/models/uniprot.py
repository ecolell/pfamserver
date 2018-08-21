from __future__ import unicode_literals
from pfamserver.database import db, Base


class Uniprot(Base):
    uniprot_acc = db.Column(
        db.UnicodeText,
        primary_key=True)
    uniprot_id = db.Column(
        db.UnicodeText,
        index=True,
        unique=True)
    description = db.Column(db.UnicodeText)
