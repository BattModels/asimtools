#!/usr/bin/env bash

# These examples explicitly specify a choice of calc_input.yaml
# If you want to use the globally set one instead, it is not necessary to 
# specify it.

# Example 1:
# This example runs multiple asimmodules one after the other, stopping on failure
# You can fix the failed job in its workdir or set the output.yaml status to 
# "discard" or "complete" and rerun asim-execute to proceed
asim-execute image_array_sim_input.yaml -c ../calc_input.yaml -e ../env_input.yaml

asim-execute image_array_key_sequence_sim_input.yaml -c ../calc_input.yaml -e ../env_input.yaml

