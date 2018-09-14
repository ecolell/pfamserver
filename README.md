# pfamserver

A python service that expose an api to the pfam database.

[![build status](https://wichi.no-ip.org/leloir/pfamserver/badges/dev/pipeline.svg)](https://wichi.no-ip.org/leloir/pfamserver/commits/dev)
[![coverage report](https://wichi.no-ip.org/leloir/pfamserver/badges/dev/coverage.svg)](https://wichi.no-ip.org/leloir/pfamserver/commits/dev)
[![Stories in Ready](https://badge.waffle.io/ecolell/pfamserver.png?label=ready&title=Ready)](https://waffle.io/ecolell/pfamserver)

## Requirements

If you want to deploy this service on any GNU/Linux or OSX system you need to have a *mysql* database named *pfamserver* (with user *pfamserver* and password *password*), a *pfamservertest* database in the same postgres engine (with the same user as owner) and the python package *virtualenvwrapper* installed.

    $ sudo mysql -u root
    > GRANT ALL PRIVILEGES ON *.* TO 'root'@'%' IDENTIFIED BY 'password';
    > GRANT ALL ON *.* TO 'root'@'%';
    > \q

First, if you are in ubuntu you could install `virtualenvwrapper` through:

    pip install --user virtualenvwrapper

And append a the next 3 lines at the bottom of the `.bashrc` or `.zshrc` file;

Configure

    export WORKON_HOME=$HOME/.virtualenvs
    export PROJECT_HOME=$HOME/Devel
    source /usr/local/bin/virtualenvwrapper.sh

Then, you should download [this repository](https://wichi.no-ip.org/leloir/pfamserver) and execute:

    mkvirtualenv pfamserver
    workon pfamserver
    pip install --upgrade pip
    pip install --upgrade distribute
    pip install -r requirements_development.txt

If you are deploying to production you should use:

    pip install -r requirements.txt

On Ubuntu Desktop there are some other libraries not installed by default which may need to be installed to run this server. Use the next command to automate the installation of the additional C libraries:

    sudo apt-get install zlibc curl libssl0.9.8 libbz2-dev libxslt-dev libxml-dev

To install the hmmer library and be able to recover a pfam alignment in stockholm format you should execute:

    flask library hmmer install -v 3.2.1
    flask library hmmer stockholm_index -v 31.0

To install the pfamscan library to be able to get a list of pfam accessions from a protein sequence you should execute:

    flask library pfamscan install
    flask library pfamscan index -v 31.0

To deploy the database structure you should execute:

    flask db shrinked download -v 31.0
    flask db shrinked install -v 31.0
    flask db data build_cache -v 31.0 

On the frontend side, it was tested with:

    node@6.11.0
    npm@5.1.0
    yarn@0.27.5

To get those versions on osx, you need to run:

    sudo port install node npm3

And to get those version on a ubuntu server, you need to run:

    sudo apt-get install nodejs npm

Then in both systems you need to run:

    sudo npm install -g n
    sudo n 6.11.0
    sudo npm install -g npm@5.1.0
    sudo npm install -g yarn@0.27.5
    npm install .

## Shrinking process

To download a specific database version from the [ebi.ac.uk] server, use the commands:

    workon pfamserver
    pip install -r requirements.txt
    flask db structure build -v 31.0
    flask db data download -v 31.0
    flask db data load -v 31.0

To shrink it and pack it (it could take a day):

    flask db data shrink -v 31.0
    flask db dump -v 31.0
    flask db pack_dump -v 31.0

To make the shrink available:

    1. Upload the file into the proper `Google Drive account`.
    2. Update the `pfamserver/commands/db.py` file. Set the reationship between `pfam_version` and `google_drive_object_id` into the `versions` dictionary.
    3. Push the changes to the main repository.
    4. Test it using:

        flask db shrinked download -v 31.0

## Testing

To test the project you should use the command:

    FLASK_ENV=testing py.test

NOTE: Remember to set SQLALCHEMY_DATABASE_URI to the right version of the pfam database.

If you want to help us or report an issue join to us through our [issue tracker](https://wichi.no-ip.org/leloir/pfamserver/issues).

## Running

To run the project you should use the command:

    FLASK_ENV=developing honcho start web webpack

NOTE: Remember to set SQLALCHEMY_DATABASE_URI to the right version of the pfam database.

## Example

A json API is available through [https://localhost:5000/api/v0](http://localhost:5000/api/v0).

There are some **special queries** that join multiple tables and accept some specific parameters. This queries respond with a json dictionary with the "query" and "output" keys.

First, it is posible to obtain a list of **pfamA** registers from a **uniprot_id** (eg. [http://localhost:5000/api/query/pfam_uniprot/egfr_human](http://localhost:5000/api/query/pfam_uniprot/egfr_human)).

Also, it is available a way to recover a list of **PDB** registers from a "**sequence description**" (eg. [http://localhost:5000/api/query/pdb_sequencedescription/egfr_human,57,168](http://localhost:5000/api/query/pdb_sequencedescription/egfr_human,57,168)) where the 3 parameters are a **uniprot_id**, a **seq_start** and **seq_end** inside that uniprot register.

Last, it is posible to obtain an MSA in stockholm format through a pfamA_acc or pfamA_id (eg. [http://localhost:5000/api/query/stockholm_pfam/piwi](http://localhost:5000/api/query/stockholm_pfam/piwi)). The **output** value is *zip compressed* and then *base 64 encoded*, to optimize transport over the network.

Please, check [http://localhost:5000/](http://localhost:5000/) for more updated examples.

## About

This software is developed by [LELOIR](http://leloir.org.ar/). You can contact us to [eloy.colell@gmail.com](mailto:eloy.colell@gmail.com).
