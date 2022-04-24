FROM alpine/make:4.2.1 as make
FROM python:3.8.10-slim AS backend_base
LABEL maintainer="Eloy Adonis Colell <eloy.colell@gmail.com>"

ENV PYTHONUNBUFFERED True
COPY --from=make /usr/bin/make /usr/bin/make

RUN mkdir -p /home/pfamserver/stage/

WORKDIR /home/pfamserver/stage

RUN pip install --upgrade pip==21.1.3
COPY requirements.txt ./
RUN pip install --no-cache-dir -r requirements.txt

FROM backend_base as backend_web

CMD exec gunicorn --bind :$PORT --workers 1 --threads 8 --timeout 30 --keep-alive 30 --error-logfile=- --access-logfile=- --worker-class=gevent wsgi:app
FROM backend_web as backend_dev

COPY requirements_dev.txt ./
RUN pip install --no-cache-dir -r requirements_dev.txt
