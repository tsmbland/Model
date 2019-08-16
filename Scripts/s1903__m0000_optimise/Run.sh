#!/usr/bin/env bash

#SBATCH --array=0-10000

module load Python/3.5.2-foss-2016b

srun python Run.py $SLURM_ARRAY_TASK_ID