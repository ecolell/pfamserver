from .extensions import db
import sqlalchemy as sa


class TimestampMixin(object):
    """Mixin to add timestamp columns to SQLAlchemy models"""

    created = sa.Column(sa.types.DateTime(timezone=True), default=sa.func.now())
    updated = sa.Column(sa.types.DateTime(timezone=True), default=sa.func.now(), onupdate=sa.func.now())


class Base(TimestampMixin, db.Model):
    """Base Class for SQLAlchemy models"""

    __abstract__ = True
