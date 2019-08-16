#!/usr/bin/env bash

IDa=$(sbatch --time=01:00:00 --parsable --wrap='srun python Init.py')
IDb=$(sbatch --time=01:00:00 --parsable --dependency=afterok:$IDa _run.sh)

for gen in {1..20}
do
    IDa=$(sbatch --time=01:00:00 --parsable --dependency=afterok:$IDb --wrap='srun python Init.py')
    IDb=$(sbatch --time=01:00:00 --parsable --dependency=afterok:$IDa _run.sh)
done
