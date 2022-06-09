import factory
from faker import Factory
from pfamserver.database import db
from pfamserver.models import (
    Pdb,
    PdbPfamAReg,
    PfamA,
    PfamARegFullSignificant,
    Pfamseq,
    Uniprot,
    UniprotRegFull,
)

faker = Factory.create()


class UniprotFactory(factory.alchemy.SQLAlchemyModelFactory):
    class Meta:
        model = Uniprot
        sqlalchemy_session = db.session


class UniprotRegFullFactory(factory.alchemy.SQLAlchemyModelFactory):
    class Meta:
        model = UniprotRegFull
        sqlalchemy_session = db.session


class PfamAFactory(factory.alchemy.SQLAlchemyModelFactory):
    class Meta:
        model = PfamA
        sqlalchemy_session = db.session


class PfamARegFullSignificantFactory(factory.alchemy.SQLAlchemyModelFactory):
    class Meta:
        model = PfamARegFullSignificant
        sqlalchemy_session = db.session


class PfamseqFactory(factory.alchemy.SQLAlchemyModelFactory):
    class Meta:
        model = Pfamseq
        sqlalchemy_session = db.session


class PdbFactory(factory.alchemy.SQLAlchemyModelFactory):
    class Meta:
        model = Pdb
        sqlalchemy_session = db.session


class PdbPfamARegFactory(factory.alchemy.SQLAlchemyModelFactory):
    class Meta:
        model = PdbPfamAReg
        sqlalchemy_session = db.session
