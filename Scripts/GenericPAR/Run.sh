#!/usr/bin/env bash

#SBATCH --array=4-4
#SBATCH --time=50:00:00
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=32

srun python kant.py $SLURM_ARRAY_TASK_ID