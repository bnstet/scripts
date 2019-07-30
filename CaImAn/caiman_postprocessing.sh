#!/bin/bash
#SBATCH --partition=cpu_short
#SBATCH --job-name=caiman_postprocess
#SBATCH --mem=4GB
#SBATCH --time=0-00:05:00
#SBATCH --tasks=1
#SBATCH --cpus-per-task=2
#SBATCH --nodes=1
#SBATCH -o /gpfs/home/stetlb01/logs/postproc_%j.log
#SBATCH -e /gpfs/home/stetlb01/scripts/CaImAn/logs/postproc_%j.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user=stetlb01@nyulangone.org


export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1


python caiman_postprocessing.py "$@"




