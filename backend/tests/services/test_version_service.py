from pfamserver.services import version_service as service


def test_version(app, current_version):
    with app.app_context():
        db_name = app.config.get("SQLALCHEMY_DATABASE_URI").split("/")[-1]
        assert db_name[:8] == current_version.replace(".", "_")
        assert service.version() == current_version
