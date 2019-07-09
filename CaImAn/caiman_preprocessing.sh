#!/bin/bash
#SBATCH --partition=cpu_short
#SBATCH --job-name=caiman_preprocess
#SBATCH --mem=50GB
#SBATCH --time=0-01:00:00
#SBATCH --tasks=1
#SBATCH --cpus-per-task=4
#SBATCH --nodes=1
#SBATCH -o ~/scripts/CaImAn/logs/preproc_%j.log
#SBATCH -e ~/scripts/CaImAn/logs/preproc_%j.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user=stetlb01@nyulangone.org

# input text file
infile=$1

echo ""
echo "********************"
echo "Running CaImAn preprocessing on video $invid with output $outfile"
echo "********************"
echo ""

export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1


finallog=logs/final_preproc_log_"$SLURM_JOB_ID".txt
printf "$(date)  Task $SLURM_JOB_ID starting analysis of videos in $infile \n" >> $finallog


python caiman_preprocessing.py $infile
rc=$?
if [[ $rc != 0 ]]; then printf "$(date)  Task $SLURM_JOB_ID failed with exit code $rc \n" >> $finallog; exit $rc; fi
printf "$(date)  Task $SLURM_JOB_ID completed \n" >> $finallog

exit


