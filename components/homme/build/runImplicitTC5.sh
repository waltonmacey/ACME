#!/bin/bash

#SBATCH --nodes=4
#SBATCH --time=24:00:00
#SBATCH --account=fy150039
#SBATCH --partition=ec
#SBATCH --mail-type=ALL
#SBATCH --job-name=hommeSWTC5

export Trilinos_DIR=/ascldap/users/pabosle/trilinos-intel-kokkosSerial/install

export NETCDFDIR=$HOME/netcdf-4.3.2-intel
export PNETCDFDIR=$HOME/pnetcdf-1.6.1
export HDF5DIR=$HOME/hdf5-1.8.12-intel

export HOMME_ROOT=/ascldap/users/pabosle/ACME/components/homme

export HOMME_EXE=$HOMME_ROOT/build/src/swim/swim 

export RUN_DIR=/fscratch/pabosle/hommeImplicit

cp $HOMME_ROOT/build/trilinosOptions.xml $RUN_DIR/.

cd $RUN_DIR
#mkdir movies

rm -rf movies/swtc5*.nc

date
mpirun -np 32 $HOMME_EXE < $HOMME_ROOT/build/swtc5-implicit.nl | tee swimOutput.txt
date

