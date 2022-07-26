# EMBL-EBI released a new PFAM version!!

Since **pfamserver** only use only 7 tables from PFAM database, the process consist of building structure (1), downloading data files (2), loading data files into mysql database (3), cropping not used columns from database (4), construct join table to improve performance (5), and dump the destiled database into a compact file to simplify distribution into the many servers around there (6).

The used tables are:

    * pfamA
    * pfamseq
    * uniprot
    * pdb
    * pdb_pfamA_reg
    * uniprot_reg_full
    * pfamA_reg_full_significant

Since it is a good practise to keep raw downloaded files from previous versions, each version has a folde r inside backend folder.


## 1. Preparing the database structure

The first step is to prepare the database structure. To do so, a software developer should use the **testingnet** database available on the docker-compose.dev.yml file.

The command

    make db-structure

should download an **\<table\>.sql.gz** file, extract it, apply **\<table\>.sql** structure to the database and remove it; for each table in the list.


## 2. Downloading the database data

The second step is to download the PFAM database content from the FTP server. Since some of the tables contains big amounts of data:

    * pfamA (9 MB)
    * pfamseq (18 GB)
    * uniprot (57 GB)
    * pdb (8 MB)
    * pdb_pfamA_reg (7 MB)
    * uniprot_reg_full (6 GB)
    * pfamA_reg_full_significant (4 GB)

it is required to have at least 85 GB of free disk space to be able to download, and around 400 GB more in order to uncompress each file and inject it into the mysql database.

On the other hand, take into account that the download could take a few days since with good internet connections the download speed could be limited by the FTP server.

The command

    make db-data-download


should download an **\<table\>.txt.gz** file for each table in the list.


## 3. Loading the downloaded content into mysql

The third step is to uncompress, load into the mysql table and remove uncompressed files. Since uniprot.txt uncompressed take (160 GB on disk and even more into the Mysql DB), you probably needs an extra 430GB in order to be able to load the 7 tables into the mysql database.

The command

    make db-data-load


should ontain a **\<table\>.txt** from a **\<table\>.txt.gz**, load it into the right table, and remove the uncompressed to optimize disk space.


## 4. Cropping not used columns

The forth step is to remove not used columns from mysql database.

The command

    make db-data-cropp

leave only useful columns and reduce database size making the useful part of PFAM more portable.


## 5. Construct join table to improve performance (cache table)

The fifth step is to create a view to improve the performance when try to get the list of sequences (with and without PDB) for a specific pfam.

The command

    make db-data-cache


## 6. Dump destiled data, compact and upload into Google Drive.

The last step is to extract the shrinked database, compact it and upload into Mistic2's Google Drive, and keep the link into this repository.

The command

    make db-data-pack

generates the **./db/pfam<version>.sql.bz2** file that should be uploaded into Google Drive.

Last, Makefile should be updated using the new GOOGLE_DRIVE_ID and PFAM_VERSION.
