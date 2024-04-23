#!/usr/bin/env bash

# These examples explicitly specify a choice of calc_input.yaml
# If you want to use the globally set one instead, it is not necessary to 
# specify it.

# Example 1: This example relaxes the cell to a hydrostatic external pressure
# of 0.05eV/Ang^3
asim-execute symmetric_cell_relax_sim_input.yaml -c ../calc_input.yaml -e ../env_input.yaml

