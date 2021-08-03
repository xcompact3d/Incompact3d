#!/bin/bash
#
# Default partition is compute
#
# Name of the job
#SBATCH --job-name=channel
#
# Number of node(s)
#SBATCH --nodes=1
#
# Number of tasks per node
#SBATCH --tasks-per-node=32
#
# Runtime of this jobs is 6 hours
#SBATCH --time=6:00:00
#
# Send SIGUSR2 240 seconds before the end
#SBATCH --signal=SIGUSR2@240

# Clear the environment from any previously loaded modules
module purge > /dev/null 2>&1

# Load the module environment suitable for the job
module load compiler/intel/2020.4.304 mpi/intel/2020.4.304

# Intel is looking for Slurm PMI
export I_MPI_PMI_LIBRARY=/usr/lib64/libpmi.so

# run the job
srun -n${SLURM_NTASKS} ./xcompact3d
