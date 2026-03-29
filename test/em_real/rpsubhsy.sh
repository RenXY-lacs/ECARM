#!/bin/bash
#SBATCH -p amd_256
#SBATCH -n 1
#SBATCH -n 1
source /public1/soft/modules/module.sh
module load wrf/4.2
srun ./wrf-pp.exe >& xwrf1-pp.log
touch WRF_OK
