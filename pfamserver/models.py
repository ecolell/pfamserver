# from application import app
from database import engine, scoped_db as db
from sqlalchemy.ext.automap import automap_base
# from sqlalchemy import MetaData
from autoupdate import loader
from sqlalchemy import Table, Column, String, MetaData, Integer, Boolean, literal_column
from sqlalchemy import and_, types
from sqlalchemy.orm import mapper
from sqlalchemy.ext.declarative import declared_attr
from sqlalchemy.ext.automap import generate_relationship, interfaces
from sqlalchemy.sql.functions import concat
from sqlalchemy.sql.expression import cast
import re
import inflect


class Base(object):

    @declared_attr
    def __tablename__(cls):
        return cls.__name__

    __table_args__ = {'mysql_engine': 'MyISAM'}


Base = automap_base(cls=Base)
Base.query = db.query_property()


def camelize_classname(base, tablename, table):
    "Produce a 'camelized' class name, e.g. "
    "'words_and_underscores' -> 'WordsAndUnderscores'"
    return str(tablename[0].upper() +
               re.sub(r'_([a-z])', lambda m: m.group(1).upper(),
                      tablename[1:]))


def pluralize_collection(base, local_cls, referred_cls, constraint):
    "Produce an 'uncamelized', 'pluralized' class name, e.g. "
    "'SomeTerm' -> 'some_terms'"
    referred_name = referred_cls.__name__
    uncamelized = re.sub(r'[A-Z]',
                         lambda m: "_%s" % m.group(0).lower(),
                         referred_name)[1:]
    pluralized = _pluralizer.plural(uncamelized)
    return pluralized


def gen_relationship(base, direction, return_fn,
                     attrname, local_cls, referred_cls, **kw):
    if direction is interfaces.ONETOMANY:
        kw['cascade'] = 'all, delete-orphan'
        kw['passive_deletes'] = True
        # make use of the built-in function to actually
        # return
        # the result.
        return generate_relationship(base, direction, return_fn,
                                     attrname, local_cls, referred_cls, **kw)


_pluralizer = inflect.engine()
Base.prepare(engine, reflect=True,
             classname_for_table=camelize_classname,
             name_for_collection_relationship=pluralize_collection,
             generate_relationship=gen_relationship)

classes = dict(Base.classes.items())
globals().update(classes)
classes = classes.values()
print "--> Coverage: ", len(classes), "/", len(loader.tables)


class Pfamjoinpfamseq(object):
    pass

def create_aux_pfamA_pfamseqid():
    from sqlalchemy.orm import aliased

    tablename = "pfamjoinpfamseq"
    m = MetaData(engine)
    if not engine.dialect.has_table(engine, tablename):
        # t.drop(engine) # to delete/drop table
        t = Table(tablename, m, Column('pfamseq_id', String(40), primary_key=True),
                                Column('pfamA_acc', String(7), primary_key=True, index=True),
                                Column('pfamseq_acc', String(40)),
                                Column('has_pdb', Boolean, unique=False, default=False))
        t.create(engine)

        pfams = db.query(PfamA.pfamA_acc).all()
        for pf in pfams:
            print pf[0]
            pfam = pf[0]

            query = db.query(concat(Pfamseq.pfamseq_id, '/',
                                           cast(PfamARegFullSignificant.seq_start, types.Unicode), '-',
                                           cast(PfamARegFullSignificant.seq_end, types.Unicode)),
                                    PfamARegFullSignificant.pfamA_acc,
                                    PfamARegFullSignificant.pfamseq_acc)

            query = query.join(PfamARegFullSignificant,
                            and_(PfamARegFullSignificant.in_full == 1,
                                 Pfamseq.pfamseq_acc == PfamARegFullSignificant.pfamseq_acc,
                                 PfamARegFullSignificant.pfamA_acc == pfam))

            ### add column with pdb
            subquery2 = db.query(PdbPfamAReg.pfamseq_acc)
            subquery2 = subquery2.filter(PdbPfamAReg.pfamA_acc == pfam).distinct().subquery()
            query_pdb = query.filter(PfamARegFullSignificant.pfamseq_acc == subquery2.c.pfamseq_acc)
            subquery_pdb = query_pdb.subquery()

            query = query.filter(PfamARegFullSignificant.pfamseq_acc.notin_(subquery2))

            # print query.add_columns(literal_column("0").label("has_pdb")).distinct().all()
            # print "################"
            # print "################"
            # print query_pdb.add_columns(literal_column("1").label("has_pdb")).distinct().all()

            query = query.add_columns(literal_column("0").label("has_pdb")).distinct()
            query_pdb = query_pdb.add_columns(literal_column("1").label("has_pdb")).distinct()
            query_union = query.union(query_pdb)
            engine.execute(
                t.insert().values(tuple(query_union.all())))
        return ["Succesfully created pfamjoinpfamseq"]

    else:
        pfamjoinpfamseq = Table(tablename, m, autoload=True )
        mapper(Pfamjoinpfamseq, pfamjoinpfamseq)
        return ["pfamjoinpfamseq table exists, creating mapper"]

_output = create_aux_pfamA_pfamseqid()
print(_output)
