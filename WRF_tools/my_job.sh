#!/bin/bash

#SBATCH -A "WRF"
#SBATCH -J WRF
#SBATCH -p high
#SBATCH --ntasks=90
#SBATCH --exclusive

#Environment
module load netcdf/intel/4.4.4
module load intel/mpi/64/2020/2.54

export FI_PROVIDER=tcp

#Run command
mpirun -np 90 ./wrf.exe

