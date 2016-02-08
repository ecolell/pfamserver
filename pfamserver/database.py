from application import app
from sqlalchemy import create_engine
from sqlalchemy.orm import scoped_session, sessionmaker

engine = create_engine(app.config['SQLALCHEMY_DATABASE_URI'])
Session = sessionmaker(autocommit=False, autoflush=False, bind=engine)
Session.configure(bind=engine)
scoped_db = scoped_session(Session)
app.config['SESSION_SQLALCHEMY'] = scoped_db()
