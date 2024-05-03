#!/usr/bin/env bash

# These examples explicitly specify a choice of calc_input.yaml
# If you want to use the globally set one instead, it is not necessary to 
# specify it.

# Example 1:
# This example runs the same simulation for Cu using different calc_id values
# See sim_array examples for an alternative way to do this
asim-execute calc_array_calc_id_sim_input.yaml -c ../calc_input.yaml -e ../env_input.yaml

# Example 2:
# This example runs the same elastic constant simulation for Cu using different calc_id values
asim-execute calc_array_elastic_constant_sim_input.yaml -c ../calc_input.yaml -e ../env_input.yaml

# Example 3:
# This example runs the same simulation on Cu using different Lennard-Jones sigma parameters
asim-execute calc_array_sigma_sim_input.yaml -c ../calc_input.yaml -e ../env_input.yaml

# Example 4:
# This example runs the same simulation on Cu using different Lennard-Jones sigma parameters
# while simultaneously looping over a second epsilon argument
asim-execute calc_array_secondary_sim_input.yaml -c ../calc_input.yaml -e ../env_input.yaml
