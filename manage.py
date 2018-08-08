#!/usr/bin/env python

import os
from pfamserver import create_app
from pfamserver.database import db
from flask_script import Manager, Server, Shell


application = create_app()
app = application


def make_shell_context():
    return {"app": app, "db": db}


manager = Manager(app)
manager.add_command('runserver', Server(host='0.0.0.0', use_reloader=True,
                                        use_debugger=app.config['DEBUG']))
manager.add_command('shell', Shell(make_context=make_shell_context))


if __name__ == '__main__':
    manager.run()
