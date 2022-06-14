#!make
# include .env
# export $(shell sed 's/=.*//' .env)

assert-command-present = $(if $(shell which $1),,$(error '$1' missing and needed for this build))

# SHELL += -x
PFAM_VERSION=35.0
HMMER_VERSION=3.2.1

MAKE=make
# UID:=$(shell id -u):$(shell id -g)
# Hardcode uid and gid to avoid gitlab ci set it with root id (or something weird).
UID:=1000:1000
DC_BASE:=CURRENT_UID=$(UID) docker-compose -f "docker-compose.yml"
DC:=$(DC_BASE) --project-name=$(PROJECT_NAME)
DC_DEV:=CURRENT_UID=$(UID) docker-compose -f "docker-compose.dev.yml"
DOCKER:=docker run --rm -v "$(PWD)/backend:/work"
DBASH:=$(DOCKER) bash:4.4

# Dev initialization
dev-init:
	@echo "You need to manually install:"
	@echo "1. docker-compose"
	@echo "2. Set a virtual environment with:"
	@echo "   mkvirtualenv pfamserver --python=`which python3.8`"
	@echo "3. Setup pre-commit."
	@echo "   curl https://pre-commit.com/install-local.py | python -"
	@echo "   pre-commit install"


# Docker images
hmmer:
	$(DBASH) rm /work/Pfam$(PFAM_VERSION)/hmmer-$(HMMER_VERSION).tar.gz
	$(DBASH) wget http://eddylab.org/software/hmmer/hmmer-$(HMMER_VERSION).tar.gz -P /work/Pfam$(PFAM_VERSION)
	$(DBASH) tar xvzf /work/Pfam$(PFAM_VERSION)/hmmer-$(HMMER_VERSION).tar.gz -C /work/Pfam$(PFAM_VERSION)
	$(DOCKER) -w /work/Pfam$(PFAM_VERSION)/hmmer-$(HMMER_VERSION) \
		gcc:4.9 ./configure
	$(DOCKER) -w /work/Pfam$(PFAM_VERSION)/hmmer-$(HMMER_VERSION) \
		gcc:4.9 make


preflight-stockholm: hmmer
	$(DBASH) wget -c ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam$(PFAM_VERSION)/Pfam-A.hmm.gz -P /work/Pfam$(PFAM_VERSION)
	$(DBASH) wget -c ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam$(PFAM_VERSION)/Pfam-A.hmm.dat.gz -P /work/Pfam$(PFAM_VERSION)
	$(DBASH) gunzip -fk /work/Pfam$(PFAM_VERSION)/Pfam-A.hmm.gz
	$(DBASH) gunzip -fk /work/Pfam$(PFAM_VERSION)/Pfam-A.hmm.dat.gz
	$(DBASH) rm -f /work/Pfam$(PFAM_VERSION)/Pfam-A.hmm.h3i
	$(DBASH) rm -f /work/Pfam$(PFAM_VERSION)/Pfam-A.hmm.h3m
	$(DOCKER) gcc:4.9 /work/Pfam$(PFAM_VERSION)/hmmer-$(HMMER_VERSION)/src/hmmpress /work/Pfam$(PFAM_VERSION)/Pfam-A.hmm
	$(DBASH) wget -c ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam$(PFAM_VERSION)/Pfam-A.full.gz -P /work/Pfam$(PFAM_VERSION)
	$(DBASH) gunzip -fk /work/Pfam$(PFAM_VERSION)/Pfam-A.full.gz
	$(DBASH) sed -i -E "s/(#=GF AC   [A-Z0-9]+)\\.(.+)/\\1\\n#=GF DC   Revision: \\2/g" /work/Pfam$(PFAM_VERSION)/Pfam-A.full


hmm-index: #preflight-stockholm
	$(DBASH) rm -f /work/Pfam$(PFAM_VERSION)/Pfam-A.full.ssi
	$(DOCKER) \
		gcc:4.9 /work/Pfam$(PFAM_VERSION)/hmmer-$(HMMER_VERSION)/easel/miniapps/esl-afetch --index /work/Pfam$(PFAM_VERSION)/Pfam-A.full


pfamscan: hmm-index
	$(DBASH) mkdir -p /work/Pfam$(PFAM_VERSION)/PfamScan
	$(DBASH) wget -c ftp://ftp.ebi.ac.uk/pub/databases/Pfam/Tools/PfamScan.tar.gz -P /work/Pfam$(PFAM_VERSION)
	$(DBASH) tar xvzf /work/Pfam$(PFAM_VERSION)/PfamScan.tar.gz -C /work/Pfam$(PFAM_VERSION)
	$(DBASH) mkdir -p /work/tmp
	$(DOCKER) -w /work/Pfam$(PFAM_VERSION)/ \
		bash:4.4 ln -s -f ./hmmer-$(HMMER_VERSION)/src/hmmscan hmmscan


pre-flight: pfamscan


docker-build-dev:
	$(DC_DEV) --verbose --log-level DEBUG build web


docker-build: # pre-flight
	# $(MAKE) -C backend extract-requirements
	$(DC) build web


docker-upload:
	-docker pull ecolell/pfamserver-web:latest
	-docker pull ecolell/pfamserver-web:latest-dev
	docker build -t "ecolell/pfamserver-web:latest" --cache-from "ecolell/pfamserver-web:latest" --target backend_base backend
	docker build -t ecolell/pfamserver-web:latest-dev --cache-from ecolell/pfamserver-web:latest-dev --target backend_dev backend
	docker login --username=ecolell
	docker push ecolell/pfamserver-web:latest
	docker push ecolell/pfamserver-web:latest-dev

# Deployment targets

pipeline-database-test:
	FLASK_APP=backend/flasky.py MIGRATION_DIR=backend/pfamserver/models/migrations $(MAKE) -C backend pipeline-database-test

# Testing targets

pipeline-backend-test: # docker-build-dev
	mkdir -p db/mysql_test backend/tmp
	$(DC_DEV) up -d db
	$(DC_DEV) run -u root -w "/home/pfamserver/stage" -e FLASK_APP=/home/pfamserver/stage -e FLASK_ENV=testing web py.test -s | tee pytest-coverage.txt
	$(DC_DEV) down

pipeline-backend-mypy:
	$(DC_DEV) run --rm -w "/home/pfamserver/stage" -e FLASK_APP=/home/pfamserver/stage web mypy pfamserver tests

pipeline-backend-security:
	$(DC_DEV) run --rm -w "/home/pfamserver/stage" -e FLASK_APP=/home/pfamserver/stage web bandit -x pfamserver/command -r pfamserver

pipeline-backend-safety:
	$(DC_DEV) run --rm -w "/home/pfamserver/stage" -e FLASK_APP=/home/pfamserver/stage web safety check --full-report

pipeline-backend-quality:
	$(DC_DEV) run --rm -w "/home/pfamserver/stage" -e FLASK_APP=/home/pfamserver/stage web pydocstyle pfamserver
	$(DC_DEV) run --rm -w "/home/pfamserver/stage" -e FLASK_APP=/home/pfamserver/stage web radon cc -nb -a pfamserver

pipeline-backend: pipeline-backend-mypy pipeline-backend-test pipeline-backend-safety pipeline-backend-security pipeline-backend-quality

pipeline-test-e2e: release-dev maintenance-off
	$(DC) run --rm e2e

# Development

stop:
	$(DC) stop

clean:
	$(DC) down -v --remove-orphans

# Release commands

rel-start:
	@echo -n "--> Releasing version "
	@cat VERSION.txt

rel-end:
	@echo "--> Release ready"

setup-libraries:
	$(DC) run -e PFAM_VERSION=$(PFAM_VERSION) -e HMMER_VERSION=$(HMMER_VERSION) --entrypoint="make" --rm web setup-libraries


shrinked-download:
	$(DC) run -e PFAM_VERSION=$(PFAM_VERSION) --entrypoint="make" --rm web-dev shrinked-download

shrinked-install:
	$(DC) run -e PFAM_VERSION=$(PFAM_VERSION) --entrypoint="make" --use-aliases --rm web-dev shrinked-install

shrinked-build-cache:
	$(DC) run -e PFAM_VERSION=$(PFAM_VERSION) --entrypoint="make" --use-aliases --rm web-dev shrinked-build-cache

setup-shrinked: shinked-download shrinked-install shrinked-build-cache

up: setup-libraries
	$(DC) up -d db
	@sleep 3;

	$(DC) up -d web nginx
	# $(DC) scale web=4

up-db-admin:
	$(DC) up -d phpmyadmin

down-db-admin:
	$(DC) stop phpmyadmin
	$(DC) rm phpmyadmin

bootup: rel-start clean docker-build up setup-shrinked rel-end
blue-green-release:
	./deploy.sh

up-traefik:
	CURRENT_UID=$(UID) docker-compose -f "docker-compose.traefik.yml" up -d

down-traefik:
	CURRENT_UID=$(UID) docker-compose -f "docker-compose.traefik.yml" down

ps:
	@(CURRENT_UID=$(UID) docker-compose -f "docker-compose.traefik.yml" ps)
	@echo "\n"
	@($(DC_BASE) --project-name=green ps)
	@echo "\n"
	@($(DC_BASE) --project-name=blue ps)

tlogs:
	@(CURRENT_UID=$(UID) docker-compose -f "docker-compose.traefik.yml" logs --tail=3 -f)

# Cleanup and utils

logs:
	$(DC) logs -f nginx web

cleanup-docker:
	docker stop $(shell docker ps -a -q)
	docker kill $(shell docker ps -q)
	docker rm $(shell docker ps -a -q)
	docker rmi $(shell docker images -q)

docker-cleanup-all:
	docker system prune -a

check-commands:
	$(call assert-command-present, docker)
	$(call assert-command-present, docker-compose)

check-env:
	env
$(DC_DEV) up -d -w "/home/pfamserver/stage" -e FLASK_APP=/home/pfamserver/stage -e FLASK_ENV=development db
