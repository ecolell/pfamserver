from application import app
from sqlalchemy import create_engine
from sqlalchemy.orm import scoped_session, sessionmaker
from autoupdate.core import Manager


version = Manager().get_versions()[0]
engine = create_engine(app.config['SQLALCHEMY_DATABASE_URI'] + version,
                       pool_size=20, max_overflow=100)
Session = sessionmaker(autocommit=False, autoflush=False, bind=engine)
Session.configure(bind=engine)
scoped_db = scoped_session(Session)
app.config['SESSION_SQLALCHEMY'] = scoped_db()


@app.teardown_appcontext
def shutdown_session(exception=None):
        scoped_db.remove()
