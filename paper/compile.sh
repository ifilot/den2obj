#!/bin/bash

# Perform Docker compilation under Windows using Git BASH

# provide the path to the location where the paper is stored; note that under
# Windows Git Bash, two slashes need to be used instead of one

if [ -d "/c/PROGRAMMING/CPP/den2obj/paper" ]
then
    PAPERPATH="//c//PROGRAMMING//CPP//den2obj//paper"
elif [ -d "/d/PROGRAMMING/CPP/den2obj/paper" ]
then
    PAPERPATH="//d//PROGRAMMING//CPP//den2obj//paper"
else
    echo "Cannot find compilation path"
    exit -1
fi

docker run --rm \
--volume $PAPERPATH:/data \
--user $(id -u):$(id -g) \
--env JOURNAL=joss \
openjournals/inara
