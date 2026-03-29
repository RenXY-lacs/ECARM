#!/bin/bash
#SBATCH -p amd_256
#SBATCH -N 1
#SBATCH -n 64
#SBATCH -J rxynewwrf
source /public23/soft/modules/module.sh
module load wrf/4.2
srun ./wrf.exe >& xwrf1.log
touch WRF_OK
