#!/usr/bin/env bash

# These examples explicitly specify a choice of env_input.yaml and calc_input.yaml
# If you want to use the globally set ones instead, it is not necessary to specify
# them.

# Run the chain of asimmodules in a mix of interactive and inline jobs
asim-execute chained_srun_sim_input.yaml -c ../calc_input.yaml -e ../env_input.yaml

# Run the chain of asimmodule modules in batch jobs with dependencies
asim-execute chained_batch_sim_input.yaml -c ../calc_input.yaml -e ../env_input.yaml

# Run a chain that internally submits an array, connected by job depndencies
asim-execute eos_chained_array_sim_input.yaml -c ../calc_input.yaml -e ../env_input.yaml
