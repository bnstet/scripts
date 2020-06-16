#!/bin/bash
#SBATCH -c 4
#SBATCH --mem=25G
#SBATCH -p cpu_short
#SBATCH -t 00:30:00
#SBATCH -o logs/gen_log_%A_%a.log
#SBATCH -e logs/gen_log_%A_%a.log


bash -c "$@"
