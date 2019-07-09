#!/bin/bash
#SBATCH --nodes=1
#SBATCH --partition=cpu_short
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --time=00:20:00
#SBATCH --mem=120GB
#SBATCH --job-name=cell_series_analysis
#SBATCH --mail-type=END
#SBATCH -o /gpfs/home/stetlb01/logs/cell_series_%j.log 

module purge
module load matlab/R2018a

loadSettings=$1
fpathPat=$2
h5Pat=$3
stackPat=$4
saveDataPat=$5
figPat=$6
cmPat=$7

{
  runCommand="addpath('/gpfs/home/stetlb01/Holography_Analysis');addpath('/gpfs/home/stetlb01/JG_Functions');CellSeriesAnalysis($loadSettings,'$fpathPat','$h5Pat','$stackPat','$saveDataPat','$figPat','$cmPat');exit"
  echo $runCommand
  matlab -nodisplay -r "$runCommand"
     }  2>&1

exit
