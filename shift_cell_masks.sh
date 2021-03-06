#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --time=0:20:00
#SBATCH --mem=20GB
#SBATCH --job-name=mask_shift
#SBATCH --mail-type=END
#SBATCH -o /gpfs/home/stetlb01/logs/shift_%A_%a.log
#SBATCH -e /gpfs/home/stetlb01/logs/shift_%A_%a.log

module purge
module load matlab/R2018a

## sourcef (1st arg): source template
## targetf (2nd arg): target template
## sourcemaskdir (3rd arg) source mask directory
## targetmaskdir (4th arg) output mask directory


sourcef=$1
targetf=$2
sourcemaskdir=$3
targetmaskdir=$4


export SLURM_DIR=/gpfs/scratch/stetlb01/${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}
mkdir -p $SLURM_DIR


echo "Matlab command:  " "addpath('/gpfs/home/stetlb01/JG_Functions');addpath('/gpfs/home/stetlb01/normcorre-matlab');shift_masks('$sourcef','$targetf','$sourcemaskdir','$targetmaskdir');exit"

{
    matlab -nodisplay -r "addpath('/gpfs/home/stetlb01/JG_Functions');addpath('/gpfs/home/stetlb01/normcorre-matlab');shift_masks('$sourcef','$targetf','$sourcemaskdir','$targetmaskdir');exit"
}

rm -r $SLURM_DIR

exit

