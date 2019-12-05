#!/bin/bash

case $# in
    3)
        sudo docker run --rm -v $1:$1 -v $2:$2 -v $2/_FEATURE/PSORTB:/tmp/results e4ong1031/vaxign-ml:latest python3.6 VaxignML.py -i $1 -o $2 -t $3
        ;;
    4)
        sudo docker run --rm -v $1:$1 -v $2:$2 -v $2/_FEATURE/PSORTB:/tmp/results e4ong1031/vaxign-ml:latest python3.6 VaxignML.py -i $1 -o $2 -t $3 -s $4
        ;;
esac