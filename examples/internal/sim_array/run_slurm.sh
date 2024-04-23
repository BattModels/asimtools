#!/usr/bin/env bash

# These examples explicitly specify a choice of calc_input.yaml
# If you want to use the globally set ones instead, it is not necessary to 
# specify them.

# Example 1:
# This example runs the same calculator on FCC Cu with different
# lattice parameters (a). 

asim-execute sim_array_batch_lattice_parameters_sim_input.yaml -c ../calc_input.yaml -e ../env_input.yaml

# Example 2:
# This example runs the same calculator on Cu with different
# crystal structures
asim-execute sim_array_batch_crystalstructure_sim_input.yaml -c ../calc_input.yaml -e ../env_input.yaml

# Example 3:
# This example runs the same Cu FCC structure with different calculators.
# Notice that you can individually specify the env_id for each calculator,
# this allows for example to have one be DFT and another be a ML potential
# which can have very different hardware and environment requirements.
# For an alternative way to do this that allows iterating over the internal
# calc_input parameters, see the workflows.calc_array asimmodule
asim-execute sim_array_batch_calc_id_sim_input.yaml -c ../calc_input.yaml -e ../env_input.yaml

# Example 4:
# This example runs the same calculator on different input structure files.
# It also autmatically names the result directories (ids) corresponding to the
# name of the file. You can do this with any input file, not just structures
asim-execute sim_array_batch_image_file_sim_input.yaml -c ../calc_input.yaml -e ../env_input.yaml
