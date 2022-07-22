#!/bin/bash

#SBATCH --nodes=1

# Choose a high-memory node
#SBATCH --constraint=skylake

# Set the numbre of tasks per node (=MPI tasks)
#SBATCH --ntasks-per-node=1
# Set number of threads per-rank (=OpenMP threads)
#SBATCH --cpus-per-task=40

# Wall clock time limit
#SBATCH --time=00:30:00

# Load dependencies
module load matlab/R2020a

# Ensure OMP threads match CPUs
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# Ensure LD library path can find compiled libraries at runtime
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/mnt/home/skailasa/lib

# Run driver program
matlab -nodisplay -nosplash -nodesktop < ie_fmm3dbie_sound_hard.m > ie_fmm3dbie_sound_hard.output
