#!/usr/bin/env bash

# This runs all the examples that use slurm and will launch a number of very
# short jobs, if any fail, feel free to check issues on github or
# submit a new one

source clean_all.sh

RUNFILE="run.sh"

for d in */; do
    cd $d
    if test -f ${RUNFILE}; then
        source run_slurm.sh
    fi
    cd ../
done
