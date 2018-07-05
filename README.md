# pfamserver

A python service that expose an api to the pfam database.

[![build status](https://wichi.no-ip.org/leloir/pfamserver/badges/master/build.svg)](https://wichi.no-ip.org/leloir/pfamserver/commits/master)
[![coverage report](https://wichi.no-ip.org/leloir/pfamserver/badges/master/coverage.svg)](https://wichi.no-ip.org/leloir/pfamserver/commits/master)
[![Stories in Ready](https://badge.waffle.io/ecolell/pfamserver.png?label=ready&title=Ready)](https://waffle.io/ecolell/pfamserver)


Requirements
============

### old
Also, you need to have the **wget** package installed and access (username and password) to a mysql database.

Last, if you cloned the repository and you want to create a boot script on Linux you should execute:

    $ sudo make create_boot_script

### new

If you want to deploy this service on any GNU/Linux or OSX system you need to have the *mysql* engine installed (with user *admin* and password *password*), and the python package *virtualenvwrapper* installed.

First, if you are in ubuntu you could install `virtualenvwrapper` through:

    $ pip install --user virtualenvwrapper

And append a the next 3 lines at the bottom of the `.bashrc` or `.zshrc` file;

Configure

    export WORKON_HOME=$HOME/.virtualenvs
    export PROJECT_HOME=$HOME/Devel
    source /usr/local/bin/virtualenvwrapper.sh

Then, you should download [this repository](https://wichi.no-ip.org/leloir/pfamserver) and execute:

    $ mkvirtualenv pfamserver
    $ workon pfamserver
    $ pip install --upgrade pip
    $ pip install --upgrade distribute
    $ pip install -r requirements_development.txt

If you are deploying to production you should use:

    $ pip install -r requirements.txt

On Ubuntu Desktop there are some other libraries not installed by default which may need to be installed to run this server. Use the next command to automate the installation of the additional C libraries:

    $ sudo apt-get install zlibc curl libssl0.9.8 libbz2-dev libxslt-dev libxml-dev

To deploy the database structure you should execute:

    $ workon pfamserver
    $ ./manage.py db upgrade && ./manage.py seed

Then, you should install **RabbitMQ** to simplify the process on **MacOS** use the next Makefile commands:

    $ sudo port install openssl
    $ make rmq-osx-download
    $ make rmq-osx-start
    $ make rmq-osx-clean

and on a **Debian-like** GNU/Linux use these other ones:

    $ sudo apt-get install openssl
    $ make rmq-linux-download
    $ make rmq-linux-start
    $ make rmq-linux-clean

On the frontend side, it was tested with:

    node@6.11.0
    npm@5.1.0
    yarn@0.27.5

To get those versions on osx, you need to run:

    $ sudo port install node npm3

And to get those version on a ubuntu server, you need to run:

    $ sudo apt-get install nodejs npm

Then in both systems you need to run:

    $ sudo npm install -g n
    $ sudo n 6.11.0
    $ sudo npm install -g npm@5.1.0
    $ sudo npm install -g yarn@0.27.5
    $ yarn install && yarn run build


Testing
=======

To test all the project you should use the command:

    $ py.test

If you want to help us or report an issue join to us through our [issue tracker](https://wichi.no-ip.org/leloir/pfamserver/issues).


Running
=======

To run the server you need to execute:

    $ DEBUG=True DB=root:root@localhost:3306/pfamserver HOST=0.0.0.0 PORT=5001 python -m pfamserver

After start, the server it is going to check for new updates at night in background. The update is composed of the download and deployment of the last version of the file *Pfam-A.full.gz* (approximately 13GB) and download and deployment of all the files inside the *database_files* folder (approximately 77GB), so be patient it can take a while.

At midnight, when the deploy of the last update is ready (because was runned on background) it should restart the machine. You should ensure the previous command is executed after the server boots up so it should stay update automatically without any human intervention.


Example
=======

A json API is available through [https://localhost:5001/api/](http://localhost:5001/api/), or use the graphical interface to discover some pfams database tables and apply some filters (access to [https://localhost:5001/admin](http://localhost:5001/admin)).

There are some **special queries** that join multiple tables and accept some specific parameters. This queries respond with a json dictionary with the "query" and "output" keys.

First, it is posible to obtain a list of **pfamA** registers from a **uniprot_id** (eg. [http://localhost:5001/api/query/pfam_uniprot/egfr_human](http://localhost:5001/api/query/pfam_uniprot/egfr_human)).

Also, it is available a way to recover a list of **PDB** registers from a "**sequence description**" (eg. [http://localhost:5001/api/query/pdb_sequencedescription/egfr_human,57,168](http://localhost:5001/api/query/pdb_sequencedescription/egfr_human,57,168)) where the 3 parameters are a **uniprot_id**, a **seq_start** and **seq_end** inside that uniprot register.

Then, it is accesible an alternative to obtain a **small PDB image** (base 64 enconded after compressed) from an specific **pdb_id** (eg. [http://localhost:5001/api/query/pdbimage_pdb/1IVO](http://localhost:5001/api/query/pdbimage_pdb/1IVO)).

Last, it is posible to obtain an MSA in stockholm format through a pfamA_acc or pfamA_id (eg. [http://localhost:5001/api/query/stockholm_pfam/piwi](http://localhost:5001/api/query/stockholm_pfam/piwi)). The **output** value is *zip compressed* and then *base 64 encoded*, to optimize transport over the network.

Please, check [http://localhost:5001/](http://localhost:5001/) for more updated examples.


About
=====

This software is developed by [LELOIR](https://leloir.org.ar/). You can contact us to [eloy.colell@gmail.com](mailto:eloy.colell@gmail.com).
