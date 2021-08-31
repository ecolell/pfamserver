api:
	gunicorn wsgi:app -R --error-logfile=- --access-logfile=- --reload --worker-class=gevent --bind=0.0.0.0:5000

python-requirements-dev:
	echo "To update python, libraries or algorithms versions need to rebuild the docker image."

pipeline-database-test:
	FLASK_ENV=development flask db upgrade
	FLASK_ENV=development flask db downgrade

pipeline-test:
	FLASK_ENV=development py.test -s -v
	pip install mypy
	mypy .

lint-mypy:
	pip install mypy
	mypy .

pipeline-security:
	FLASK_ENV=development bandit -x pfamserver/command -r pfamserver

pipeline-safety:
	FLASK_ENV=development safety check --full-report

pipeline-quality:
	FLASK_ENV=development pycodestyle pfamserver tests
	FLASK_ENV=development pydocstyle pfamserver
	FLASK_ENV=development radon cc -nb -a pfamserver

db-migrate:
	FLASK_ENV=development flask db migrate

db-upgrade:
	FLASK_ENV=development flask db upgrade

db-downgrade:
	FLASK_ENV=development flask db downgrade

update-setup-py:
	pip install poetry2setup
	poetry2setup > setup.py
