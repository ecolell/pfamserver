from __future__ import unicode_literals
from builtins import str as text

import factory
from faker import Factory

from pfamserver.database import db
from pfamserver.models import Uniprot

faker = Factory.create()


class UniprotFactory(factory.alchemy.SQLAlchemyModelFactory):

    class Meta:
        model = Uniprot
        sqlalchemy_session = db.session