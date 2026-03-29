#!/bin/bash
#SBATCH -p amd_512
#SBATCH -n 1
#SBATCH -n 1
source /public1/soft/modules/module.sh
module load wrf/4.2
srun ./ideal.exe >& wrf_ideal.log
touch WRF_OK
