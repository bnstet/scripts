#!/bin/bash
#SBATCH --partition=cpu_short
#SBATCH --job-name=caiman_mc
#SBATCH --mem=50GB
#SBATCH --time=0-4:00:00
#SBATCH --tasks=1
#SBATCH --cpus-per-task=20
#SBATCH --nodes=1
#SBATCH -o logs/mc_%j.log
#SBATCH -e logs/mc_%j.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user=stetlb01@nyulangone.org



echo ""
echo "********************"
echo "Running CaImAn motion correct on video $1"
echo "********************"
echo ""

export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1


finallog=logs/final_log_"$SLURM_ARRAY_JOB_ID".txt
printf "$(date)  Task $SLURM_ARRAY_JOB_ID starting analysis of $vid \n" >> $finallog



python motion_correct_only.py "$@"
rc=$?
if [[ $rc != 0 ]]; then printf "$(date)  Task $SLURM_JOB_TASK_ID failed with exit code $rc \n" >> $finallog; exit $rc; fi
printf "$(date)  Task $SLURM_ARRAY_JOB_ID completed \n" >> $finallog

exit


