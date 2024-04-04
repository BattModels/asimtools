#!/usr/bin/env bash

# This runs all the examples and might take a while, if any fail, feel free to 
# check issues on github or submit a new one

source clean_all.sh

RUNFILE="run.sh"

for d in */; do
    cd $d
    if test -f ${RUNFILE}; then
        source ${RUNFILE}
    fi
    cd ../
done
