"""Database abstract classes."""
from typing import Any

import sqlalchemy as sqla

from .extensions import db

Model: Any = db.Model


class TimestampMixin:
    """Mixin to add timestamp columns to SQLAlchemy models."""

    created = sqla.Column(sqla.types.DateTime(timezone=True), default=sqla.func.now())
    updated = sqla.Column(
        sqla.types.DateTime(timezone=True),
        default=sqla.func.now(),
        onupdate=sqla.func.now(),
    )


class Base(TimestampMixin, Model):
    """Base Class for SQLAlchemy models."""

    __abstract__ = True
