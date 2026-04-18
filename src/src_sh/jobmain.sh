#!/bin/bash

#SBATCH --job-name sic_interface4 
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=48
#SBATCH --time=36:30:00
#SBATCH --account=iontransport
#SBATCH --error=std.err_%j
#SBATCH --output=std.out_%j
#SBATCH --partition medmem

# Load LAMMPS module
module purge
module use /nopt/nrel/apps/modules/centos77/modulefiles/
module load mpich
module load intel
module load lammps//062322-intel-mpich

# Change directory and signal job start
cd $SLURM_SUBMIT_DIR
echo "begin job.."
echo $PWD

# Create output directories
mkdir -p outdir

#echo "nve/limit starting ..."
#srun --mpi=pmi2 lmp -in in.init -e screen
#wait

#echo "langevin starting ..."
#srun --mpi=pmi2 lmp -in in.langevin -e screen
#wait

echo "NVT run starting ..."
srun --mpi=pmi2 lmp -in in.nvt -e screen
echo "Calculations completed ..."
