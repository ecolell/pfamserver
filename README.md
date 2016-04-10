[![Stories in Ready](https://badge.waffle.io/ecolell/pfamserver.png?label=ready&title=Ready)](https://waffle.io/ecolell/pfamserver)
pfamserver
==========

This is a PFAM server to deploy in an intranet.


Requirements
============

If you want to deploy this service on any GNU/Linux or OSX system you just need to execute:

    $ pip install pfamserver

If you want to improve this library, you should download the [github repository](https://github.com/ecolell/pfamserver) and execute:

    $ make deploy

On Ubuntu Desktop there are some other libraries not installed by default (zlibc curl libssl0.9.8 libbz2-dev libxslt-dev libxml-dev) which may need to be installed to use these library. Use the next command to automate the installation of the additional C libraries:

    $ make ubuntu deploy

Also, you need to have the **wget** package installed and access (username and password) of a mysql database.


Testing
=======

To test all the project you should use the command:

    $ make test

If you want to help us or report an issue join to us through the [Github issue tracker](https://github.com/ecolell/pfamserver/issues).


Running
=======

To run the server you need to execute:

    $ DEBUG=True DB=root:root@localhost:3306/pfamserver HOST=0.0.0.0 PORT=5001 python -m pfamserver

After start, the server it is going to check for new updates at night in background. The update is composed of the download and deployment of the last version of the file *Pfam-A.full.gz* (approximately 13GB) and download and deployment of all the files inside the *database_files* folder (approximately 77GB), so be patient it can take a while.

At midnight, when the deploy of the last update is ready (because was runned on background) it should restart the machine. You should ensure the previous command is executed after the server boots up so it should stay update automatically without any human intervention.


Example
=======

A json API is available through [http://localhost:5001/api/](http://localhost:5001/api/), or use the graphical interface to discover some pfams database tables and apply some filters (access to [http://localhost:5001/admin](http://localhost:5001/admin)).

There are some **special queries** that join multiple tables and accept some specific parameters. This queries respond with a json dictionary with the "query" and "output" keys.

First, it is posible to obtain a list of **pfamA** registers from a **uniprot_id** (eg. [http://localhost:5001/api/query/pfam_uniprot/egfr_human](http://localhost:5001/api/query/pfam_uniprot/egfr_human)).

Also, it is available a way to recover a list of **PDB** registers from a "**sequence description**" (eg. [http://localhost:5001/api/query/pdb_sequencedescription/egfr_human,57,168](http://localhost:5001/api/query/pdb_sequencedescription/egfr_human,57,168)) where the 3 parameters are a **uniprot_id**, a **seq_start** and **seq_end** inside that uniprot register.

Then, it is accesible an alternative to obtain a **small PDB image** (base 64 enconded after compressed) from an specific **pdb_id** (eg. [http://localhost:5001/api/query/pdbimage_pdb/1IVO](http://localhost:5001/api/query/pdbimage_pdb/1IVO)). 

Last, it is posible to obtain an MSA in stockholm format through a pfamA_acc or pfamA_id (eg. [http://localhost:5001/api/query/stockholm_pfam/piwi](http://localhost:5001/api/query/stockholm_pfam/piwi)). The **output** value is *zip compressed* and then *base 64 encoded*, to optimize transport over the network.

Please, check [http://localhost:5001/](http://localhost:5001/) for more updated examples.


About
=====

This software is developed by [LELOIR](http://leloir.org.ar/). You can contact us to [eloy.colell@gmail.com](mailto:eloy.colell@gmail.com).
