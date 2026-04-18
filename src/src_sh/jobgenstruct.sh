#!/bin/bash

#SBATCH --job-name genstruct_interf5
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=23:30:00
#SBATCH --account=iontransport
#SBATCH --error=std.err_%j
#SBATCH --output=std.out_%j
#SBATCH --partition shared

# Load LAMMPS module
module purge
module use /nopt/nrel/apps/modules/centos77/modulefiles/
module load mpich
module load intel

# Change directory and signal job start
cd $SLURM_SUBMIT_DIR
echo "begin job.."
echo $PWD


echo "generating initial structure ..."
~/tools/packmol/packmol < gen_layer_structure.inp
echo "Calculations completed ..."
