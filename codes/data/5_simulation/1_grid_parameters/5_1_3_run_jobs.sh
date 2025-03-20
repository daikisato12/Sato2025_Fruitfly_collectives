#!/bin/sh
#SBATCH -a 1-1100
#SBATCH --mem-per-cpu 4g
#SBATCH -J job_grid_simulation_group
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

ulimit -s unlimited
echo running on `hostname`
echo starting at
date

CONDITION=`cat ./grid_simulation_group_parameters.txt | awk -v line=$SLURM_ARRAY_TASK_ID '{if (NR == line) print $0}'`
echo ${CONDITION}

python 5_1_1_grid_simulation_group.py 100 ${CONDITION}
Rscript 5_1_0_make_dataset.R grid_simulation_group 100 ${CONDITION}

echo ending at
date