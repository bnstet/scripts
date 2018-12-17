#!/bin/bash
#SBATCH --partition=gpu4_short
#SBATCH --job-name=caiman_pipeline
#SBATCH --mem=30GB
#SBATCH --time=0-02:00:00
#SBATCH --tasks=1
#SBATCH --cpus-per-task=16
#SBATCH --nodes=1
#SBATCH --array=0
#SBATCH -o logs/caiman_%A_%j_%a.log
#SBATCH -e logs/caiman_%A_%j_%a.log
#SBATCH --mail-type=ALL,ARRAY_TASKS
#SBATCH --mail-user=stetlb01@nyulangone.org


invid=$1
outfile=$2

echo ""
echo "********************"
echo "Running CaImAn pipeline on video $invid with output $outfile"
echo "********************"
echo ""

export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1


finallog=logs/final_log_"$SLURM_ARRAY_JOB_ID".txt
printf "$(date)  Task $SLURM_ARRAY_TASK_ID starting analysis of $invid \n" >> $finallog


python pipeline.py $invid $outfile --slurmid $SLURM_ARRAY_TASK_ID
rc=$?
if [[ $rc != 0 ]]; then printf "$(date)  Task $SLURM_ARRAY_TASK_ID failed with exit code $rc \n" >> $finallog; exit $rc; fi
printf "$(date)  Task $SLURM_ARRAY_TASK_ID completed \n" >> $finallog

exit


