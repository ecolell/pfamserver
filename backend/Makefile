POETRY := docker run --rm -v $(PWD)/backend:/app ecolell/poetry:3.8-slim


extract-requirements:
	# docker build -f Dockerfile.poetry --tag ecolell/poetry:3.8-slim --build-arg PORT=5000 .
	$(POETRY) export --without-hashes -f requirements.txt --output requirements.txt
	$(POETRY) export --without-hashes --dev -f requirements.txt -o requirements_dev.txt

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
