#!/bin/bash
blue=$(docker ps -f name=blue -q)
green=$(docker ps -f name=green -q)
if test -z "$blue"
then
    ENV="blue"
    OLD="green"
else
    ENV="green"
    OLD="blue"
fi

echo "Starting "$ENV" cloud"
PROJECT_NAME=$ENV make bootup

source update_traefik_config.sh

echo "Stopping "$OLD" cloud"
PROJECT_NAME=$OLD make clean
