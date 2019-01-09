#!/bin/bash
#SBATCH --partition=cpu_short
#SBATCH --job-name=caiman_mc
#SBATCH --mem=10GB
#SBATCH --time=0-01:00:00
#SBATCH --tasks=1
#SBATCH --cpus-per-task=16
#SBATCH --nodes=1
#SBATCH --array=0-35
#SBATCH -o logs/mc_%A_%j_%a.log
#SBATCH -e logs/mc_%A_%j_%a.log
#SBATCH --mail-type=ALL,ARRAY_TASKS
#SBATCH --mail-user=stetlb01@nyulangone.org


# Make sure there are enough numbers in the "array" range above to hold all of the jobs.

# input arg should be a text file with the input file (directory) names

mapfile -t vid_files < $1


if [ ${#vid_files[@]} -lt 1 ]
then
    echo "ERROR: no video files provided in list." >&2
    exit
fi


if [ $SLURM_ARRAY_TASK_ID -ge  ${#vid_files[@]}  ]
then
    echo "Array task ID greater than needed. Exiting."
    exit
fi

vid=${vid_files[$SLURM_ARRAY_TASK_ID]}

echo ""
echo "********************"
echo "Running CaImAn motion correct on video $vid"
echo "********************"
echo ""

export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1


finallog=logs/final_log_"$SLURM_ARRAY_JOB_ID".txt
printf "$(date)  Task $SLURM_ARRAY_TASK_ID starting analysis of $vid \n" >> $finallog


python motion_correct_only.py $vid
rc=$?
if [[ $rc != 0 ]]; then printf "$(date)  Task $SLURM_ARRAY_TASK_ID failed with exit code $rc \n" >> $finallog; exit $rc; fi
printf "$(date)  Task $SLURM_ARRAY_TASK_ID completed \n" >> $finallog

exit


