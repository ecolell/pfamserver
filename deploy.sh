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

echo "Backup "$OLD" cloud"
PROJECT_NAME=$OLD make backup-db || true
echo "Starting "$ENV" cloud"
if test -z "$blue" && test -z "$green"
then
    PROJECT_NAME=$ENV make bootup
else
    PROJECT_NAME=$ENV make release-bootup
fi

echo "Change redirect from "$OLD" to "$ENV
curl -XPUT "http://localhost:8081/api/providers/rest" -d '{
  "http": {
    "routers": {
      "pfamserver_router": {
        "entryPoints": [
          "web"
        ],
        "service": "nginx",
        "rule": "Host(`localhost`) || Host(`pfamserver.leloir.org.ar`)"
      }
    },
    "services": {
      "nginx": {
        "loadBalancer": {
          "servers": [
            {
              "url": "http://'$ENV'_nginx_1:5001/"
            }
          ],
          "passHostHeader": true
        }
      }
    }
  }
}'

echo "Stopping "$OLD" cloud"
PROJECT_NAME=$OLD make clean
