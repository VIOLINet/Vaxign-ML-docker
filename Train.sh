#!/bin/bash

sudo docker run --rm -v $1:$1 -v $2:$2 -v $3:$3 -v $3/_FEATURE/PSORTB:/tmp/results e4ong1031/vaxign-ml:latest python3.6 lib/train.py -p $1 -n $2 -o $3 -t $4