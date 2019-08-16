#!/usr/bin/env bash

module load Python/3.5.2-foss-2016b

IDa=$(sbatch --parsable --wrap='srun python Init.py')
IDb=$(sbatch --parsable --dependency=afterok:$IDa _run.sh)

for gen in {1..99}
do
    IDa=$(sbatch --parsable --dependency=afterok:$IDb --wrap='srun python Init.py')
    IDb=$(sbatch --parsable --dependency=afterok:$IDa _run.sh)
done
