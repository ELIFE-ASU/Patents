#!/bin/bash

#SBATCH -N 1
#SBATCH -n 28
#SBATCH -J SureChemBL
#SBATCH -o slurm.test.out
#SBATCH -e slurm.test.err
#SBATCH -t 0-4:00:00

module purge
module load anaconda/py3
module load rclone/1.43

eval $"source activate my-rdkit-env"

python get_network_data.py
