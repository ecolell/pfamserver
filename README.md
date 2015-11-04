pfamserver
==========

This is a PFAM server to deploy in an intranet.


Requirements
============

If you want to deploy this service on any GNU/Linux or OSX system you just need to execute:

    $ pip install pfamserver

If you want to improve this library, you should download the [gitlab repository](https://github.com/ecolell/pfamserver) and execute:

    $ make deploy

On Ubuntu Desktop there are some other libraries not installed by default (zlibc curl libssl0.9.8 libbz2-dev libxslt-dev libxml-dev) which may need to be installed to use these library. Use the next command to automate the installation of the additional C libraries:

    $ make ubuntu deploy


Testing
=======

To test all the project you should use the command:

    $ make test

If you want to help us or report an issue join to us through the [GitLab issue tracker](https://github.com/ecolell/pfamserver/issues).


Example
=======


About
=====

This software is developed by [LELOIR](http://leloir.org.ar/). You can contact us to [eloy.colell@gmail.com](mailto:eloy.colell@gmail.com).
