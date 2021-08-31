from __future__ import unicode_literals
from pfamserver.database import db, Base


class Uniprot(Base):
    uniprot_acc = db.Column(
        db.String(10),
        primary_key=True)
    uniprot_id = db.Column(
        db.String(16),
        index=True,
        unique=True)
    description = db.Column(db.UnicodeText)
