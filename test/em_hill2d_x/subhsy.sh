#!/bin/bash
#SBATCH -p amd_256
#SBATCH -n 1
#SBATCH -n 1
source /public23/soft/modules/module.sh
module load wrf/4.2
srun ./wrf.exe >& xwrf1.log
touch WRF_OK
