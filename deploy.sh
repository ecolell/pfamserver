#!/bin/bash
blue=$(docker ps -f name=pfs-blue -q)
green=$(docker ps -f name=pfs-green -q)
if test -z "$blue"
then
    ENV="pfs-blue"
    OLD="pfs-green"
else
    ENV="pfs-green"
    OLD="pfs-blue"
fi

echo "Starting "$ENV" cloud"
PROJECT_NAME=$ENV make bootup

source update_traefik_config.sh

echo "Stopping "$OLD" cloud"
PROJECT_NAME=$OLD make clean
