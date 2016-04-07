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

    $ DEBUG=True DB=root:root@localhost:3306/pfamserver python -m pfamserver

After start, the server it is going to check for new updates at night in background. The update is composed of the download and deployment of the last version of the file *Pfam-A.full.gz* (approximately 13GB) and download and deployment of all the files inside the *database_files* folder (approximately 77GB), so be patient it can take a while.

At midnight, when the deploy of the last update is ready (because was runned on background) it should restart the machine. You should ensure the previous command is executed after the server boots up so it should stay update automatically without any human intervention.


Example
=======

You can make requests to [http://localhost:5001](http://localhost:5001), like:


    [http://localhost:5001/api/query/stockholm_pfam/piwi](http://localhost:5001/api/query/stockholm_pfam/piwi)


The response is a javascript dictionary with the "query" and "output" keys. The output value is zip compressed and then base 64 encoded, to optimize transport over the network.

Also there is a graphical interface without password to discover some pfams database tables at:


    [http://localhost:5001/admin](http://localhost:5001/admin)


Last, you should check the [http://localhost:5001/](http://localhost:5001/) for some examples.


About
=====

This software is developed by [LELOIR](http://leloir.org.ar/). You can contact us to [eloy.colell@gmail.com](mailto:eloy.colell@gmail.com).
