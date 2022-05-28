from contextlib import closing

import mysql.connector
import pytest
from pfamserver import create_app


def raw_query(query: str):
    mydb = mysql.connector.connect(
        host="db",
        user="root",
        passwd="root",
        database="sys",
    )
    cursor = mydb.cursor(buffered=True)
    cursor.execute(query)
    mydb.commit()
    mydb.disconnect()


@pytest.fixture(scope="session")
def session_app():
    raw_query("CREATE DATABASE IF NOT EXISTS Pfam35_0;")
    app = create_app()
    yield app
    raw_query("DROP DATABASE IF EXISTS Pfam35_0;")


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
    app = session_app
    app.db.create_all()
    yield app.db
    app.db.session.remove()
    engine = app.db.get_engine(app)
    metadata = app.db.Model.metadata

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
    app.db.session.begin(subtransactions=True)

    yield app.db

    app.db.session.rollback()
