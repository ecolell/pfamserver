#!make
# include .env
# export $(shell sed 's/=.*//' .env)

assert-command-present = $(if $(shell which $1),,$(error '$1' missing and needed for this build))

# SHELL += -x

MAKE=make
# UID:=$(shell id -u):$(shell id -g)
# Hardcode uid and gid to avoid gitlab ci set it with root id (or something weird).
UID:=1000:1000
DC_BASE:=CURRENT_UID=$(UID) docker-compose -f "docker-compose.yml"
DC:=$(DC_BASE) --project-name=$(PROJECT_NAME)
DC_DEV:=$(DC_BASE) run --rm -w "/home/pfamserver/stage" -e FLASK_APP=/home/pfamserver/stage --entrypoint="make" web-dev

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

docker-build:
	$(DC) build web web-dev

docker-upload:
	docker pull ecolell/pfamserver-web:latest
	docker pull ecolell/pfamserver-web:latest-dev
	docker build -t "ecolell/pfamserver-web:latest" --cache-from "ecolell/pfamserver-web:latest" --target backend_base backend
	docker build -t ecolell/pfamserver-web:latest-dev --cache-from ecolell/pfamserver-web:latest-dev --target backend_dev backend
	docker login --username=ecolell
	docker push ecolell/pfamserver-web:latest
	docker push ecolell/pfamserver-web:latest-dev

# Frontend

build-frontend:
	$(DC) run --rm frontend rm -rf /usr/app/backend/pfamserver/static/dist
	$(DC) run --rm frontend npm install --no-audit
	$(DC) run --rm frontend npm run build

# Deployment targets

pipeline-assets:
	npm $(NPM_PROXY) install --no-audit && npm run build
	mkdir -p backend/pfamserver/static
	cp -rv static/* backend/pfamserver/static/
	cp -rv backend/pfamserver/static public/

pipeline-database-test:
	FLASK_APP=backend/flasky.py MIGRATION_DIR=backend/pfamserver/models/migrations $(MAKE) -C backend pipeline-database-test

# Testing targets

mock-assets:
	-mkdir -p ./backend/pfamserver/static
	echo '{"assets":{},"publicPath":""}' > backend/pfamserver/static/manifest.json

mock-traefik-gate:
	-docker network create pfamserver_traefik_gate

unmock-traefik-gate:
	-docker network remove pfamserver_traefik_gate

pipeline-backend-test: mock-assets mock-traefik-gate
	touch backend/.env
	$(DC) run --rm -w "/home/pfamserver/stage" -e FLASK_APP=/home/pfamserver/stage --entrypoint="make" gitlab-web-dev pipeline-test
	$(MAKE) unmock-traefik-gate

pipeline-backend-security: mock-traefik-gate
	touch backend/.env
	$(DC) run --rm -w "/home/pfamserver/stage" -e FLASK_APP=/home/pfamserver/stage --entrypoint="make" gitlab-web-dev pipeline-security
	$(MAKE) unmock-traefik-gate

pipeline-backend-safety: mock-traefik-gate
	touch backend/.env
	$(DC) run --rm -w "/home/pfamserver/stage" -e FLASK_APP=/home/pfamserver/stage --entrypoint="make" gitlab-web-dev pipeline-safety
	$(MAKE) unmock-traefik-gate

pipeline-backend-quality: mock-traefik-gate
	touch backend/.env
	$(DC) run --rm -w "/home/pfamserver/stage" -e FLASK_APP=/home/pfamserver/stage --entrypoint="make" gitlab-web-dev pipeline-quality
	$(MAKE) unmock-traefik-gate

pipeline-backend: pipeline-backend-test pipeline-backend-safety pipeline-backend-security pipeline-backend-quality

pipeline-test-e2e: release-dev maintenance-off
	$(DC) run --rm e2e

# Development

stop:
	$(DC) stop

clean:
	$(DC) down -v --remove-orphans

# DB commands

load-db-examples:
	cp ./examples.sql ./backend/examples.sql
	$(DC) run --rm web flask data jobs load -f /home/pfamserver/stage/examples.sql

build-db: load-db-examples

backup-db: dump-db
	@echo "--> Extracting DB from docker compose."

dump-db:
	$(DC) run --rm web flask data jobs dump -f /home/pfamserver/stage/dumped.sql

# Release commands

rel-start:
	@echo -n "--> Releasing version "
	@cat VERSION.txt

rel-end:
	@echo "--> Release ready"

up: # dump-db
	$(DC) up -d db
	@sleep 3;

	$(DC) up -d web nginx
	# $(DC) scale web=4

restore-db:
	$(DC) run --rm web flask data jobs load -f /home/pfamserver/stage/dumped.sql

bootup: rel-start clean docker-build build-frontend up rel-end
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

# Services

celery-purge:
	$(DC_DEV) celery-purge

# Cleanup and utils

logs:
	$(DC) logs -f

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