#!/usr/bin/env bash

# These examples explicitly specify a choice of calc_input.yaml
# If you want to use the globally set ones instead, it is not necessary to specify
# them.

# Example 1: This example uses the custom defined script ase_eos.py to 
# calculate the eos for one strucuture

export ASIMTOOLS_ASIMMODULE_DIR="./" # Set the directory to start looking for assimmodule
asim-execute ase_eos_sim_input.yaml -c calc_input.yaml -e ../env_input.yaml

# Example 2: This recreates the ASE example for finding the eos for fcc
# metals, but is now completely parallelized and can be used for any
# structure with cubic symmetry

asim-execute array_ase_eos_sim_input.yaml -c calc_input.yaml -e ../env_input.yaml
