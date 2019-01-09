#!/bin/bash
#SBATCH --partition=gpu4_short
#SBATCH --job-name=caiman_pipeline
#SBATCH --mem=30GB
#SBATCH --time=0-06:00:00
#SBATCH --tasks=1
#SBATCH --cpus-per-task=16
#SBATCH --nodes=1
#SBATCH --array=0-10
#SBATCH -o logs/caiman_%A_%j_%a.log
#SBATCH -e logs/caiman_%A_%j_%a.log
#SBATCH --mail-type=ALL,ARRAY_TASKS
#SBATCH --mail-user=stetlb01@nyulangone.org


# Make sure there are enough numbers in the "array" range above to hold all of the jobs.


## List all video files here. 
## If the video is a .tif (single file), point directly to the .tif file. 
## If the video is a .tiff stack (many .tiff files), point to the containing folder and the algorithm will run on evert .tiff in the folder (so don't group multiple .tiff stacks in the same folder).

# Make sure there are enough numbers in the "array" range above to hold all of the jobs.

# input arg 1 should be a text file with the input file (directory) names

mapfile -t vid_files < $1

# input arg 2 should be the folder to which the outputs go

out_folder=$2



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
fname=$(basename $vid)

echo ""
echo "********************"
echo "Running CaImAn pipeline on video $vid"
echo "********************"
echo ""

export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1


finallog=logs/final_log_"$SLURM_ARRAY_JOB_ID".txt
printf "$(date)  Task $SLURM_ARRAY_TASK_ID starting analysis of $vid \n" >> $finallog


#python pipeline.py $vid results/caiman_analysis_"$SLURM_ARRAY_JOB_ID"_"$SLURM_ARRAY_TASK_ID".hdf5  --slurmid $SLURM_ARRAY_TASK_ID
python pipeline.py $vid $out_folder/caiman_$fname.hdf5
rc=$?
if [[ $rc != 0 ]]; then printf "$(date)  Task $SLURM_ARRAY_TASK_ID failed with exit code $rc \n" >> $finallog; exit $rc; fi
printf "$(date)  Task $SLURM_ARRAY_TASK_ID completed \n" >> $finallog

exit


