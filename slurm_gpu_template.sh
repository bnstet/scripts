#!/bin/bash
#SBATCH --partition=gpu4_short
#SBATCH --job-name=caiman_gpu
#SBATCH --mem=200GB
#SBATCH --time=0-01:00:00
#SBATCH --tasks=1
#SBATCH --cpus-per-task=12
#SBATCH --nodes=1
#SBATCH --gres=gpu:1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=stetlb01@nyulangone.org
#SBATCH -o "caiman_gpu.log"
#SBATCH -e "caiman_gpu.log"

srun caiman_gpu_subscript.sh
