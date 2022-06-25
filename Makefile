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
SBASH:=docker run --rm -v "$(PWD)/db:/work" bash:4.4
TABLES:=pfamA pfamseq uniprot pdb pdb_pfamA_reg uniprot_reg_full pfamA_reg_full_significant

# Dev initialization
dev-init:
	@echo "You need to manually install:"
	@echo "1. docker-compose"
	@echo "2. Set a virtual environment with:"
	@echo "   mkvirtualenv pfamserver --python=`which python3.8`"
	@echo "3. Setup pre-commit."
	@echo "   curl https://pre-commit.com/install-local.py | python -"
	@echo "   pre-commit install"


# PFAMSCAN and HMMER preparation

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


# EMBL-EBI PFAM DB destilation

MYSQL_SHELL:=docker run --rm -it -e MYSQL_ROOT_PASSWORD=root -v "/home/eloy/version/git/pfamserver/db:/work" -w /work/Pfam$(PFAM_VERSION) --network=pfamserver_testingnet mysql:8.0.26
DB_NAME:=Pfam$(subst .,_,$(PFAM_VERSION))

db-check-version:
	@$(SBASH) echo Current version: $(PFAM_VERSION)
	@$(SBASH) echo Available version: `wget -q -O - "http://ftp.ebi.ac.uk/pub/databases/Pfam/releases/?C=M;O=D" | \
			sed -n 's/.*href="Pfam\([0-9]\+\).\([0-9]\+\).*/\1.\2/p' | sort -V | tail -n 1`


db-check-size:
	@$(DC_DEV) up -d db
	@$(MYSQL_SHELL) mysql -u root -proot -h db $(DB_NAME) -e " \
	SELECT table_schema, \
        ROUND(SUM(data_length + index_length) / 1024 / 1024 / 1024, 1) 'DB Size in GB' \
        FROM information_schema.tables \
        GROUP BY table_schema HAVING table_schema='$(DB_NAME)'; \
	"

db-structure:
	@$(DC_DEV) up -d db
	@$(MYSQL_SHELL) bash -c 'echo "CREATE DATABASE IF NOT EXISTS $(DB_NAME)" | mysql -u root -proot -h db'
	@$(SBASH) mkdir -p /work/Pfam$(PFAM_VERSION)
	@for table in $(TABLES); do \
		$(SBASH) echo Building $$table structure; \
		$(SBASH) wget -c http://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam$(PFAM_VERSION)/database_files/$$table.sql.gz -O /work/Pfam$(PFAM_VERSION)/$$table.sql.gz; \
		$(SBASH) gunzip -fk /work/Pfam$(PFAM_VERSION)/$$table.sql.gz; \
		$(MYSQL_SHELL) bash -c "cat /work/Pfam$(PFAM_VERSION)/$$table.sql | mysql -u root -proot -h db $(DB_NAME)"; \
		$(SBASH) rm /work/Pfam$(PFAM_VERSION)/$$table.sql; \
	done

db-data-download:
	@for table in $(TABLES); do \
		$(SBASH) echo Downloading $$table data; \
		$(SBASH) wget -c http://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam$(PFAM_VERSION)/database_files/$$table.txt.gz -O /work/Pfam$(PFAM_VERSION)/$$table.txt.gz; \
	done


UNUSED_COLUMNS:= \
	pdb:keywords \
	pdb_pfamA_reg:auto_pdb_reg \
	pdb_pfamA_reg:pdb_start_icode \
	pdb_pfamA_reg:pdb_end_icode \
	pdb_pfamA_reg:seq_start \
	pdb_pfamA_reg:seq_end \
	pdb_pfamA_reg:hex_color \
	pfamA:previous_id \
	pfamA:author \
	pfamA:deposited_by \
	pfamA:seed_source \
	pfamA:type \
	pfamA:comment \
	pfamA:sequence_GA \
	pfamA:domain_GA \
	pfamA:sequence_TC \
	pfamA:domain_TC \
	pfamA:sequence_NC \
	pfamA:domain_NC \
	pfamA:buildMethod \
	pfamA:model_length \
	pfamA:searchMethod \
	pfamA:msv_lambda \
	pfamA:msv_mu \
	pfamA:viterbi_lambda \
	pfamA:viterbi_mu \
	pfamA:forward_lambda \
	pfamA:forward_tau \
	pfamA:num_seed \
	pfamA:version \
	pfamA:number_archs \
	pfamA:number_species \
	pfamA:number_structures \
	pfamA:number_ncbi \
	pfamA:number_meta \
	pfamA:average_length \
	pfamA:percentage_id \
	pfamA:average_coverage \
	pfamA:change_status \
	pfamA:seed_consensus \
	pfamA:full_consensus \
	pfamA:number_shuffled_hits \
	pfamA:number_uniprot \
	pfamA:rp_seed \
	pfamA:number_rp15 \
	pfamA:number_rp35 \
	pfamA:number_rp55 \
	pfamA:number_rp75 \
	pfamA_reg_full_significant:auto_pfamA_reg_full \
	pfamA_reg_full_significant:ali_start \
	pfamA_reg_full_significant:ali_end \
	pfamA_reg_full_significant:model_start \
	pfamA_reg_full_significant:model_end \
	pfamA_reg_full_significant:domain_bits_score \
	pfamA_reg_full_significant:domain_evalue_score \
	pfamA_reg_full_significant:sequence_bits_score \
	pfamA_reg_full_significant:sequence_evalue_score \
	pfamA_reg_full_significant:cigar \
	pfamA_reg_full_significant:tree_order \
	pfamA_reg_full_significant:domain_order \
	pfamseq:seq_version \
	pfamseq:crc64 \
	pfamseq:md5 \
	pfamseq:evidence \
	pfamseq:length \
	pfamseq:species \
	pfamseq:taxonomy \
	pfamseq:is_fragment \
	pfamseq:sequence \
	pfamseq:ncbi_taxid \
	pfamseq:auto_architecture \
	pfamseq:treefam_acc \
	pfamseq:swissprot \
	uniprot:seq_version \
	uniprot:crc64 \
	uniprot:md5 \
	uniprot:evidence \
	uniprot:species \
	uniprot:taxonomy \
	uniprot:is_fragment \
	uniprot:sequence \
	uniprot:ncbi_taxid \
	uniprot:ref_proteome \
	uniprot:complete_proteome \
	uniprot:treefam_acc \
	uniprot:rp15 \
	uniprot:rp35 \
	uniprot:rp55 \
	uniprot:rp75 \
	uniprot_reg_full:ali_start \
	uniprot_reg_full:ali_end \
	uniprot_reg_full:model_start \
	uniprot_reg_full:model_end \
	uniprot_reg_full:domain_bits_score \
	uniprot_reg_full:domain_evalue_score \
	uniprot_reg_full:sequence_bits_score \
	uniprot_reg_full:sequence_evalue_score


db-data-load-and-cropp:
	@$(DC_DEV) up -d db
	@$(MYSQL_SHELL) mysql -u root -proot -h db $(DB_NAME) -e"SET GLOBAL local_infile=1;"
	@for table in $(TABLES); do \
		$(SBASH) echo Loading $$table data; \
		$(SBASH) gunzip -fk /work/Pfam$(PFAM_VERSION)/$$table.txt.gz; \
		$(MYSQL_SHELL) mysql --local_infile=1 -u root -proot -h db $(DB_NAME) -e "LOAD DATA LOCAL INFILE '/work/Pfam$(PFAM_VERSION)/$$table.txt' INTO TABLE $$table CHARACTER SET latin1 COLUMNS TERMINATED BY '\\t' LINES TERMINATED BY '\\n';"; \
		$(foreach TABLECOLUMN, $(UNUSED_COLUMNS), \
			$(eval a_table = $(word 1,$(subst :, ,$(TABLECOLUMN)))) \
			$(eval column = $(word 2,$(subst :, ,$(TABLECOLUMN)))) \
			$(SBASH) echo Removing $(column) from $(a_table); \
			$(MYSQL_SHELL) mysql -u root -proot -h db $(DB_NAME) -e " \
				set @exist_Check := ( \
					select count(*) from information_schema.columns \
					where table_name='$(a_table)' \
					and '$(a_table)'='$$table' \
					and column_name='$(column)' \
					and table_schema=database() \
				) ; \
				set @sqlstmt := if(@exist_Check>0,'ALTER TABLE $(a_table) DROP COLUMN $(column);' , 'select 1') ; \
				prepare stmt from @sqlstmt ; \
				execute stmt ; \
			"; \
		) \
		$(SBASH) rm /work/Pfam$(PFAM_VERSION)/$$table.txt; \
	done
	@$(MYSQL_SHELL) mysql -u root -proot -h db $(DB_NAME) -e " \
		UPDATE uniprot SET created=updated; \
		UPDATE pfamseq SET created=updated; \
		UPDATE pfamA SET created=updated; \
	"

db-data-setcreated:
	@$(DC_DEV) up -d db
	@$(MYSQL_SHELL) mysql -u root -proot -h db $(DB_NAME) -e " \
		UPDATE uniprot SET created=updated; \
		UPDATE pfamseq SET created=updated; \
		UPDATE pfamA SET created=updated; \
	"

db-data-cropp:
	@$(DC_DEV) up -d db
	@$(foreach TABLECOLUMN, $(UNUSED_COLUMNS), \
        $(eval table = $(word 1,$(subst :, ,$(TABLECOLUMN)))) \
        $(eval column = $(word 2,$(subst :, ,$(TABLECOLUMN)))) \
		$(SBASH) echo Removing $(column) from $(table); \
		$(MYSQL_SHELL) mysql -u root -proot -h db $(DB_NAME) -e " \
			set @exist_Check := ( \
				select count(*) from information_schema.columns \
				where table_name='$(table)' \
				and column_namjkjkjjjjjkjkjje='$(column)' \
				and table_schema=database() \
			) ; \
			set @sqlstmt := if(@exist_Check>0,'ALTER TABLE $(table) DROP COLUMN $(column);' , 'select 1') ; \
			prepare stmt from @sqlstmt ; \
			execute stmt ; \
		"; \
    )


db-data-cache:
	@echo Create join table.
	@$(DC_DEV) up -d db

db-data-pack:
	@echo Dump database, gziped, and update Makefile Google Drive
	@$(DC_DEV) up -d db


# Docker images

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
