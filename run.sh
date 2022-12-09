#!/bin/bash
#SBATCH --partition=modi_short        # specify the partition to run on
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --mem=250gb                   # Job memory request
#SBATCH --time=02:00:00               # Time limit hrs:min:sec
#SBATCH --output=serial_test_%j.log   # Standard output and error log
srun singularity exec ~/modi_images/slurm-notebook-latest.sif ~/modi_mount/MutationMapper/run.sh
