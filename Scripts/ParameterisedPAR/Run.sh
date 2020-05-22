#!/usr/bin/env bash

#SBATCH --array=0-80
#SBATCH --time=50:00:00
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=32

srun python NewPF.py $SLURM_ARRAY_TASK_ID