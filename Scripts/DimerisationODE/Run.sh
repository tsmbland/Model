#!/usr/bin/env bash

#SBATCH --array=0-5
#SBATCH --time=50:00:00
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=32

srun python kon_binary.py $SLURM_ARRAY_TASK_ID