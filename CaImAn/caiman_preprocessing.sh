#!/bin/bash
#SBATCH --partition=cpu_short
#SBATCH --job-name=caiman_preprocess
#SBATCH --mem=20GB
#SBATCH --time=0-00:30:00
#SBATCH --tasks=1
#SBATCH --cpus-per-task=2
#SBATCH --nodes=1
#SBATCH -o logs/caiman_preproc_%j.log
#SBATCH -e logs/caiman_preproc_%j.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user=stetlb01@nyulangone.org


export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1


finallog=logs/final_preproc_log_"$SLURM_JOB_ID".txt
printf "$(date)  Task $SLURM_JOB_ID starting preprocessing on videos \n" >> $finallog


python caiman_preprocessing.py "$@"
rc=$?
if [[ $rc != 0 ]]; then printf "$(date)  Task $SLURM_JOB_ID failed with exit code $rc \n" >> $finallog; exit $rc; fi
printf "$(date)  Task $SLURM_JOB_ID completed \n" >> $finallog
exit

