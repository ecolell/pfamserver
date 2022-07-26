# Preparing DB for distribution

## Preparing the shrinked database

To deploy the database structure you should execute:

    flask db shrinked download -v 32.0
    flask db shrinked install -v 32.0
    flask db data build_cache -v 32.0

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
    flask db structure build -v 32.0
    flask db data download -v 32.0
    flask db data load -v 32.0

To shrink it and pack it (it could take a day):

    flask db data shrink -v 32.0
    flask db dump -v 32.0
    flask db pack_dump -v 32.0

To make the shrink available:

    1. Upload the file into the proper `Google Drive account`.
    2. Update the `pfamserver/commands/db.py` file. Set the reationship between `pfam_version` and `google_drive_object_id` into the `versions` dictionary.
    3. Push the changes to the main repository.
    4. Test it using:

        flask db shrinked download -v 32.0
