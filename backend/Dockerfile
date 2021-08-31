FROM python:3.8.10 AS backend_base
MAINTAINER Eloy Adonis Colell <eloy.colell@gmail.com>

ARG DEBIAN_FRONTEND=noninteractive

# user setup and permisssion
RUN groupadd -g 1000 pfamserver && \
    useradd -u 1000 -g pfamserver -G sudo -m -s /bin/bash pfamserver && \
    mkdir -p /home/pfamserver/.local && \
    chown -R pfamserver:pfamserver /home/pfamserver

RUN mkdir -p /home/pfamserver/requirements/
RUN (curl -L https://packagecloud.io/github/git-lfs/gpgkey | apt-key add - ) && \
    (curl -s https://packagecloud.io/install/repositories/github/git-lfs/script.deb.sh | bash) && \
    apt-get update && \
    apt-get install -y --no-install-recommends git-lfs postgresql && \
    git lfs install && \
    rm -rf /var/lib/apt/lists/*

RUN pip install --upgrade pip==21.1.3
RUN curl -sSL https://raw.githubusercontent.com/python-poetry/poetry/master/get-poetry.py -o get-poetry.py
RUN POETRY_VERSION=1.0.5 python get-poetry.py
COPY ./poetry.lock /home/pfamserver/requirements/
COPY ./pyproject.toml /home/pfamserver/requirements/

ENV PATH="${PATH}:/root/.poetry/bin"
WORKDIR /home/pfamserver/requirements
RUN poetry config virtualenvs.create false \
  && poetry install --no-dev --no-interaction --no-ansi

WORKDIR /home/pfamserver/stage

FROM backend_base as backend_web

WORKDIR /home/pfamserver/stage

ENV FLASK_ENV SQLALCHEMY_DATABASE_URI NODE_ENV SERVER_NAME FLASK_APP

EXPOSE 5000 4369

FROM backend_web as backend_dev

WORKDIR /home/pfamserver/requirements
RUN poetry install --no-interaction --no-ansi


FROM backend_worker as backend_dev_worker

WORKDIR /home/pfamserver/requirements
RUN poetry install --no-interaction --no-ansi

WORKDIR /home/pfamserver/stage

ENTRYPOINT ["make"]
CMD ["-C", "backend", "pipeline-test"]
