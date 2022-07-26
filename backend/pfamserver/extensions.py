"""Extensions builders."""

from flask_caching import Cache
from flask_sqlalchemy import SQLAlchemy
from flask_wtf.csrf import CSRFProtect

cache = Cache()
csrf = CSRFProtect()
db = SQLAlchemy()  # session_options={"autoflush": False})
