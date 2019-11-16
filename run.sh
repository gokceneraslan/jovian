#!/bin/bash

docker pull gokceneraslan/jovian
docker run -it --entrypoint '' -u jovyan -v `pwd`:/home/jovyan gokceneraslan/jovian /bin/bash run-docker.sh
