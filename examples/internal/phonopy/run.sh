#!/usr/bin/env bash

# These examples explicitly specify a choice of calc_input.yaml
# If you want to use the globally set one instead, it is not necessary to 
# specify it.

# Example 1: This chained workflow calculates the phonon dispersion and density
# of states. FOr analysis, it is best to load the phonopy_save.yaml (specified
# by the user) and use phonopy's great interface tools
# asim-execute phonon_bands_and_dos_sim_input.yaml -c ../calc_input.yaml

# Example 2: This does a full quasiharmonic approximation, predicting a number
# of thermal properties
asim-execute phonon_bands_and_dos_sim_input.yaml -c ../calc_input.yaml -e ../env_input.yaml



