#!/bin/sh -l
#SBATCH -c 32
#SBATCH -N 1
#SBATCH -q defq
#SBATCH --exclude=gpu01
python $PATH_TO_JOB
