#!/usr/bin/env bash

#SBATCH --array=0-200

srun python Run.py $SLURM_ARRAY_TASK_ID