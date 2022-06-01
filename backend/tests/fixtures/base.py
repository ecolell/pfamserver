import os
from contextlib import closing

import pytest
from pfamserver import create_app
from sqlalchemy_utils import create_database, database_exists, drop_database


@pytest.fixture(scope="session")
def session_app():
    db_uri = os.getenv(
        "SQLALCHEMY_DATABASE_URI", "mysql+pymysql://root:root@db:3306/Pfam35_0"
    )
    if not database_exists(db_uri):
        create_database(db_uri)
    print("Using {db_uri} for tests".format(db_uri=db_uri))
    app = create_app()
    yield app
    drop_database(db_uri)


@pytest.fixture(scope="class")
def app(request, session_app):
    """Establish an application context before running the tests."""
    app = session_app

    _ctx = app.test_request_context()
    _ctx.push()

    if request.cls:
        request.cls.app = app
        request.cls.client = app.test_client(use_cookies=True)

    yield app

    _ctx.pop()


@pytest.fixture(scope="class")
def require_pdb(app):
    tmp = app.config
    app.config["REQUIRE_PDB"] = True
    return tmp


@pytest.fixture(scope="class")
def no_require_pdb(app):
    tmp = app.config
    app.config["REQUIRE_PDB"] = False
    return tmp


@pytest.fixture(scope="function")
def client(app):
    return app.test_client()


@pytest.fixture(scope="class")
def testdb(session_app):
    """Establish an application context before running the tests."""
    db = session_app.db
    db.create_all()
    yield db
    db.session.remove()
    engine = db.get_engine(session_app)
    metadata = db.Model.metadata

    with closing(engine.connect()) as con:
        trans = con.begin()
        for table in reversed(metadata.sorted_tables):
            con.execute(table.delete())
        trans.commit()

    engine.dispose()


@pytest.fixture(scope="function")
def db(app, testdb):
    """
    SetUp before each test is run: push a context and use subtransactions.
    """
    testdb.session.begin(subtransactions=True)

    yield testdb

    testdb.session.rollback()
