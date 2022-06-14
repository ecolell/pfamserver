FROM perl:5.30.0 AS perl_base
FROM python:3.8.10 AS backend_base
LABEL maintainer="Eloy Adonis Colell <eloy.colell@gmail.com>"

ARG DEBIAN_FRONTEND=noninteractive
RUN apt update && apt install -y gcc
COPY --from=perl_base /usr/local/lib/perl5 /usr/local/lib/perl5
COPY --from=perl_base /usr/local/bin/perl /usr/local/bin/perl
COPY --from=perl_base /usr/local/bin/cpan /usr/local/bin/cpan
COPY --from=perl_base /usr/local/bin/cpanm /usr/local/bin/cpanm
COPY --from=perl_base /usr/bin/make /usr/bin/make
# COPY --from=perl_base /usr/bin/perl /usr/bin/perl
# COPY --from=perl_base /usr/bin/perl /usr/bin/cpanm
ENV PERL5LIB="/usr/local/lib/site_perl:$PERL5LIB"
RUN echo $PATH
RUN echo $PERL5LIB
RUN cpanm -v Moose IPC::Run

ENV PYTHONUNBUFFERED True
RUN mkdir -p /home/pfamserver/stage/

WORKDIR /home/pfamserver/stage

RUN pip install --upgrade pip==21.1.3
COPY requirements.txt ./
RUN pip install --no-cache-dir -r requirements.txt

FROM backend_base as backend_web

CMD exec gunicorn --bind :$PORT --workers 1 --threads 8 --timeout 30 --keep-alive 30 --error-logfile=- --access-logfile=- --worker-class=gevent wsgi:app
FROM backend_base as backend_dev

COPY requirements_dev.txt ./
RUN pip install --no-cache-dir -r requirements_dev.txt
