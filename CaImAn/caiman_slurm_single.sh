#!/bin/bash
#SBATCH --partition=fn_medium
#SBATCH --job-name=caiman_pipeline
#SBATCH --mem=500GB
#SBATCH --time=0-10:00:00
#SBATCH --tasks=1
#SBATCH --cpus-per-task=16
#SBATCH --nodes=1
#SBATCH -o logs/caiman_%j.log
#SBATCH -e logs/caiman_%j.log
#SBATCH --mail-type=ALL
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


finallog=logs/final_log_"$SLURM_JOB_ID".txt
printf "$(date)  Task $SLURM_JOB_ID starting analysis of $invid \n" >> $finallog


export SLURM_SUBMIT_DIR=/gpfs/scratch/stetlb01/${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}
mkdir -p $SLURM_SUBMIT_DIR



python pipeline.py "$@" --slurmid $SLURM_JOB_ID
rc=$?
if [[ $rc != 0 ]]; then printf "$(date)  Task $SLURM_JOB_ID failed with exit code $rc \n" >> $finallog; exit $rc; fi
printf "$(date)  Task $SLURM_JOB_ID completed \n" >> $finallog

rm -r $SLURM_SUBMIT_DIR

exit


