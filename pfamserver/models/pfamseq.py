from __future__ import unicode_literals
from pfamserver.database import db, Base


class Pfamseq(Base):
    pfamseq_acc = db.Column(
        db.UnicodeText,
        primary_key=True)
    pfamseq_id = db.Column(
        db.UnicodeText,
        index=True)
