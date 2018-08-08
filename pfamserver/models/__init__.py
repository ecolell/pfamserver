#from .pdb_pfam_a_reg import PdbPfamAReg


"""
# from application import app
from database import engine, scoped_db as db
from sqlalchemy.ext.automap import automap_base
# from sqlalchemy import MetaData
from autoupdate import loader
from sqlalchemy.ext.declarative import declared_attr
from sqlalchemy.ext.automap import generate_relationship, interfaces
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
"""