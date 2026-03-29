#!/bin/bash
#SBATCH -p amd_256
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -J hsybaknc 
source /public23/soft/modules/module.sh
module load wrf/4.2
srun ./wrf.exe
touch WRF_OK
