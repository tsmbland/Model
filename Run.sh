#!/usr/bin/env bash

#SBATCH --time=0
#SBATCH --cpus-per-task=32
#SBATCH --array=0-1

module load Python/3.5.2-foss-2016b

srun python -m Scripts.s180725a $SLURM_ARRAY_TASK_ID $SLURM_ARRAY_TASK_MAX