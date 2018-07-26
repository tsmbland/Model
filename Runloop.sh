#!/usr/bin/env bash


for i in {0..10}
do
    sbatch --dependency=singleton Run.sh
done