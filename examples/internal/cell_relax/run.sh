#!/usr/bin/env bash

# These examples explicitly specify a choice of calc_input.yaml
# If you want to use the globally set ones instead, it is not necessary to specify
# them.

# Example 1: This relaxes all cell degrees of freedom using ase.filter.strainfilter
# If you want to maintain symmetry, use 
# geometry_optimization.symmetric_cell_relax asimmodule
asim-execute cell_relax_sim_input.yaml -c ../calc_input.yaml -e ../env_input.yaml
