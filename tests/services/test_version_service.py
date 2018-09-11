from __future__ import unicode_literals
from pfamserver.services import version_service as service


def test_version(app):
    with app.app_context():
        assert app.config.get('SQLALCHEMY_DATABASE_URI') == 'mysql+pymysql://root:root@192.168.0.105:3306/Pfam31_0'
        assert service.version() == 'Pfam31.0'
