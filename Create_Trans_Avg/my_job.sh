#!/bin/bash

#SBATCH -A "CREATE_2021"
#SBATCH -J CREATE_2021
#SBATCH -p high
#SBATCH --ntasks=1
#SBATCH --nodelist=hpccomp102

#Environment
module load netcdf/intel/4.4.4
module load intel/mpi/64/2020/2.54

#Run command
./run_create
