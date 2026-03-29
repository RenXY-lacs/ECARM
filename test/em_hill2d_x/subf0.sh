#!/bin/bash
#SBATCH -p amd_512
#SBATCH -n 1
#SBATCH -n 1
source /public1/soft/modules/module.sh
module load wrf/4.2
srun ./wrf-hill.exe >& xwrf1.log
touch WRF_OK
