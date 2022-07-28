# pfamserver

A python service that expose an api to the pfam database.

[![build status](https://github.com/ecolell/pfamserver/actions/workflows/backend.yml/badge.svg?branch=dev)](https://github.com/ecolell/pfamserver/commits/dev)


## Example

A json API is available through [https://localhost:5000/api/v0](http://localhost:5000/api/v0).

There are some **special queries** that join multiple tables and accept some specific parameters. This queries respond with a json dictionary with the "query" and "output" keys.

First, it is posible to obtain a list of **pfamA** registers from a **uniprot_id** (eg. [http://localhost:5000/api/query/pfam_uniprot/egfr_human](http://localhost:5000/api/query/pfam_uniprot/egfr_human)).

Also, it is available a way to recover a list of **PDB** registers from a "**sequence description**" (eg. [http://localhost:5000/api/query/pdb_sequencedescription/egfr_human,57,168](http://localhost:5000/api/query/pdb_sequencedescription/egfr_human,57,168)) where the 3 parameters are a **uniprot_id**, a **seq_start** and **seq_end** inside that uniprot register.

Last, it is posible to obtain an MSA in stockholm format through a pfamA_acc or pfamA_id (eg. [http://localhost:5000/api/query/stockholm_pfam/piwi](http://localhost:5000/api/query/stockholm_pfam/piwi)). The **output** value is *zip compressed* and then *base 64 encoded*, to optimize transport over the network.

Please, check [http://localhost:5000/](http://localhost:5000/) for more updated examples.

## About

This software is developed by [LELOIR](http://leloir.org.ar/). You can contact us to [eloy.colell@gmail.com](mailto:eloy.colell@gmail.com).
