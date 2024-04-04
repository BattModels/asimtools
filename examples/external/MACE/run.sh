#!/usr/bin/env bash

# These examples explicitly specify a choice of env_input.yaml and calc_input.yaml
# If you want to use the globally set ones instead, it is not necessary to specify
# them.

asim-execute atom_relax_sim_input.yaml -c calc_input.yaml

asim-execute mace-off_atom_relax_sim_input.yaml -c calc_input.yaml
