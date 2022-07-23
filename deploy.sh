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
