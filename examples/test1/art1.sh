#!/bin/bash
#SBATCH --mem-per-cpu 4G
#SBATCH -n 8
#SBATCH --qos main
#SBATCH -J 10_11_2000atoms
#SBATCH -C intel
#SBATCH -t 01:20:00

srun --mpi=pmi2 ./bart.sh
