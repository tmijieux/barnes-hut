#!/usr/bin/env bash

#SBATCH --job-name=mijieux0
#SBATCH --output=out.0
#SBATCH --error=err.0
#SBATCH -p mistral
#SBATCH --time=02:00:00
#SBATCH --exclusive
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task 20

WORKDIR=${WORKDIR:-${HOME}/barnes-hut/cpp_bh/}
# MPIEXEC=mpiexec

cd ${WORKDIR}
. ./.module.load

./particle 8000000

