from contextlib import closing

import pytest
from pfamserver import create_app
from pfamserver.database import db as _db
from sqlalchemy import MetaData


@pytest.fixture(scope="session")
def session_app(request):
    app = create_app()
    return app


@pytest.fixture(scope="session")
def session_db(request, session_app):
    """Session-wide test database."""
    app = session_app
    _db.app = app

    def teardown():
        meta = MetaData()
        with closing(_db.engine.connect()) as con:
            trans = con.begin()
            for table in reversed(meta.sorted_tables):
                con.execute(table.delete())
            trans.commit()

    # teardown()
    _db.create_all()

    # request.addfinalizer(teardown)
    return _db


@pytest.fixture(scope="class")
def app(request, session_app):
    """Establish an application context before running the tests."""
    app = session_app

    _ctx = app.test_request_context()
    _ctx.push()

    def teardown():

        _ctx.pop()

    if request.cls:
        request.cls.app = app
        request.cls.client = app.test_client(use_cookies=True)
    # request.addfinalizer(teardown)
    return app


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
def testdb(session_db, request, session_app):
    """Establish an application context before running the tests."""
    app = session_app

    def teardown():
        app.db.session.remove()
        engine = app.db.get_engine(app)
        metadata = app.db.Model.metadata

        with closing(engine.connect()) as con:
            trans = con.begin()
            for table in reversed(metadata.sorted_tables):
                con.execute(table.delete())
            trans.commit()

        engine.dispose()

    # request.addfinalizer(teardown)
    return app


@pytest.fixture(scope="function")
def db(testdb, app, request):
    """
    SetUp before each test is run: push a context and use subtransactions.
    """

    app.db.session.begin(subtransactions=True)

    def teardown():
        """
        TearDown after each test has run: rollback any dirty changes,
        close session, pop context for cleanup.
        """
        app.db.session.rollback()
        app.db.session.close()

    request.addfinalizer(teardown)
    return app.db
