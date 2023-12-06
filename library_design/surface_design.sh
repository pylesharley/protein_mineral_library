#!/bin/bash
#SBATCH -p short
#SBATCH --mem=4g
#SBATCH -o out

CMD=$(head -n $SLURM_ARRAY_TASK_ID surface_design.list | tail -1)
exec ${CMD}

#sbatch -a 1-$(cat surface_design.list|wc -l) surface_design.sh
