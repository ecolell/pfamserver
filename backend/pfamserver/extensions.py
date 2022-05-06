"""Extensions builders."""

from flask_caching import Cache
from flask_sqlalchemy import SQLAlchemy
from flask_wtf.csrf import CSRFProtect

cache = Cache(config={"CACHE_TYPE": "simple"})
csrf = CSRFProtect()
db = SQLAlchemy()  # session_options={"autoflush": False})
