#!/bin/bash
#
#SBATCH --mail-type=ALL
#SBATCH --cpus-per-task=5
#SBATCH -t 00:05:00

module load MATLAB
matlab -nodisplay -nosplash < topR.m

export MCR_CACHE_ROOT=$TMPDIR