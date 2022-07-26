"""Mistic2 Click commands."""
from pfamserver.extensions import db
from pfamserver import create_app, register_cli

app = create_app()


@app.shell_context_processor
def make_shell_context():
    """Make shell context."""
    return dict(app=app, db=db)
