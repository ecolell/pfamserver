from __future__ import unicode_literals
from pfamserver.services import version_service as service


def test_version(app, current_version):
    with app.app_context():
        assert app.config.get('SQLALCHEMY_DATABASE_URI')[-8:] == current_version.replace('.', '_')
        assert service.version() == current_version
