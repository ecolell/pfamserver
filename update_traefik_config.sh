#!/bin/bash
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
