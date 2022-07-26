# Working Dockerized

If you want to work with this service on any GNU/Linux or OSX system you need to have **docker compose** support and be able to run Makefile's targets.

## Build

The build stage, take care of downloading specified versions of **HMMER**, **Pfam-A.hmm** and **Pfam-A.full** files. It also build the backend image and prepare those files to be used by the API. Running

    $ make docker-build


## Test

After build stage, a minimal test on the code could be runned using:

    $ make pipeline-backend

 This is not going to use those external files nor extended data from Pfam database, but it is going to allow a fast check over the sourcecode features.


## Release

Last, we could deploy the service running

    $ make up-traefik
    $ make blue-green-release

while **traefik** expose the port 5001 to be used by the service. The **blue-green-deploy** mount a continue release system to reduce downtimes.


## Support

If you want to help us or report an issue join to us through our [issue tracker](https://wichi.no-ip.org/leloir/pfamserver/issues).
