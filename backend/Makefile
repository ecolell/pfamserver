api:
	gunicorn wsgi:app -R --error-logfile=- --access-logfile=- --reload --worker-class=gevent --bind=0.0.0.0:5000

python-requirements-dev:
	echo "To update python, libraries or algorithms versions need to rebuild the docker image."

pipeline-test:
	FLASK_ENV=development py.test -s -v
	pip install mypy
	mypy .

lint-mypy:
	pip install mypy
	mypy pfamserver tests

pipeline-security:
	FLASK_ENV=development bandit -x pfamserver/command -r pfamserver

pipeline-safety:
	FLASK_ENV=development safety check --full-report

pipeline-quality:
	FLASK_ENV=development pycodestyle pfamserver tests
	FLASK_ENV=development pydocstyle pfamserver
	FLASK_ENV=development radon cc -nb -a pfamserver

setup-libraries:
	FLASK_ENV=development flask library hmmer install -v ${HMMER_VERSION}
	FLASK_ENV=development flask library hmmer stockholm_index -v ${PFAM_VERSION}
	FLASK_ENV=development flask library pfamscan install
	FLASK_ENV=development flask library pfamscan index -v ${PFAM_VERSION}

shrinked-download:
	FLASK_ENV=development flask db shrinked download -v ${PFAM_VERSION}

shrinked-install:
	FLASK_ENV=development flask db shrinked install -v ${PFAM_VERSION}

shrinked-build-cache:
	FLASK_ENV=development flask db data build-cache -v ${PFAM_VERSION}

update-setup-py:
	pip install poetry2setup
	poetry2setup > setup.py
