#!/usr/bin/env bash

# This runs all the examples and might take a while, if any fail, feel free to 
# check issues on github or submit a new one

RUNFILE="run.sh"

for d in */; do
    cd $d
    rm -r *results
    cd ../
done
