#!/bin/bash
#SBATCH -p amd_512
#SBATCH -n 1
#SBATCH -n 1
source /public1/soft/modules/module.sh
module load wrf/4.2
srun ./compile em_squall2d_y >& xCPL2.log
touch WRFcpl_OK
