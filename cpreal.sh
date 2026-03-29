#!/bin/bash
#SBATCH -p v6_384
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -J cpreal
source /public5/soft/modules/module.sh
module load wrf/4.2
srun ./compile em_real -j 6  >& xCPL2.log
touch WRFcpl_OK
