#!/bin/bash
#SBATCH -c 4
#SBATCH --mem=25G
#SBATCH -p cpu_short
#SBATCH -t 00:30:00


bash -c "$@"
