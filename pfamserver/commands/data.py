

# verify old process to create the required primary keys and undefined indexes to speed up the queries.
# sudo mysql -u root -p pfamserver -e "LOAD DATA LOCAL INFILE '/home/eloy/version/git/pfamserver/dumps/pfamA.txt' INTO TABLE pfamA CHARACTER SET latin1 COLUMNS TERMINATED BY '\t' LINES TERMINATED BY '\n';"