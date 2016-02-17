[![Beerpay](http://test.beerpay.io/ecolell/pfamserver/badge.svg?style=flat-square)](http://test.beerpay.io/ecolell/pfamserver)
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

Also, you need to have the **wget** package installed.


Testing
=======

To test all the project you should use the command:

    $ make test

If you want to help us or report an issue join to us through the [Github issue tracker](https://github.com/ecolell/pfamserver/issues).


Running
=======

To run the server you need to execute:

    $ python -m pfamserver

Before the start, the server it is going to update (when it is outdated) to the last version of the file Pfam-A.full.gz (approximately 13GB), so be patient it can take a while.


Example
=======

You can make requests to [http://localhost:5001](http://localhost:5001), like:


    [http://localhost:5001/api/query/piwi](http://localhost:5001/api/query/piwi)


About
=====

This software is developed by [LELOIR](http://leloir.org.ar/). You can contact us to [eloy.colell@gmail.com](mailto:eloy.colell@gmail.com).
