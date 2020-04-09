#!/usr/bin/env bash

#SBATCH --array=1-12
#SBATCH --time=50:00:00
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=32

srun python Run.py $SLURM_ARRAY_TASK_ID