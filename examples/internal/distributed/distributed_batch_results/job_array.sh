#!/usr/bin/sh
#SBATCH -N 1
#SBATCH -n 2
#SBATCH -p cpu
if [[ ! -z ${SLURM_ARRAY_TASK_ID} ]]; then
    fls=( id-* )
    WORKDIR=${fls[${SLURM_ARRAY_TASK_ID}]}
fi

CUR_DIR=`pwd`
cd ${WORKDIR}

source ~/.bashrc
conda activate asimtools

echo "LAUNCHDIR: ${CUR_DIR}"
echo "WORKDIR: ${WORKDIR}"
echo "Job started on `hostname` at `date`"
  asim-run sim_input.yaml --calc=calc_input.yaml   

cd ${CUR_DIR}
echo "Job ended at `date`"