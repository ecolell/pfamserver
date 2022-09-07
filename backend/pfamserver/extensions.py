"""Extensions builders."""

from flask_caching import Cache
from flask_sqlalchemy import SQLAlchemy
from flask_wtf.csrf import CSRFProtect
from flask import request


cache = Cache()
csrf = CSRFProtect()
db = SQLAlchemy()  # session_options={"autoflush": False})


def make_cache_key(*args, **kwargs):
    """Logic to build redis-cache's key."""
    path = request.path
    args = str(list(request.args.items()))
    return path + args
