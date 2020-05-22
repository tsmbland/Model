#!/usr/bin/env bash

#SBATCH --array=0-35
#SBATCH --time=50:00:00
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=32

srun python x_kon0_linear.py $SLURM_ARRAY_TASK_ID