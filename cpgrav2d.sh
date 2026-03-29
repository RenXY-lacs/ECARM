#!/bin/bash
#SBATCH -p v6_384
#SBATCH -n 1
#SBATCH -n 1
source /public1/soft/modules/module.sh
module load wrf/4.2
srun ./compile em_grav2d_x >& xCPL-grav2d.log
touch WRFcpl_OK
