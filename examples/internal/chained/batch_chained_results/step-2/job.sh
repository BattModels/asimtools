#!/usr/bin/sh

#SBATCH -N 1
#SBATCH -n 2
#SBATCH -p cpu
#SBATCH -J step-2

source ~/.bashrc
conda activate asimtools

echo "Job started on `hostname` at `date`"
  asim-run sim_input.yaml --calc=calc_input.yaml   
echo "Job ended at `date`"