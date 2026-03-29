#!/bin/bash
#SBATCH -p amd_256
#SBATCH -N 1
#SBATCH -n 8
source /public23/soft/modules/module.sh
module load wrf/4.2
srun ./real.exe >& xreal1.log
