from application import app
from sqlalchemy import create_engine
from sqlalchemy.orm import scoped_session, sessionmaker

engine = create_engine(app.config['SQLALCHEMY_DATABASE_URI'], convert_unicode=True)
db = sessionmaker(autocommit=False, autoflush=False, bind=engine)
scoped_db = scoped_session(db)
app.config['SESSION_SQLALCHEMY'] = db
