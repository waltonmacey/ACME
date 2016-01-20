#!/bin/bash

export Trilinos_DIR=/ascldap/users/pabosle/trilinos-intel-kokkosSerial/install

export NETCDFDIR=$HOME/netcdf-4.3.2-intel
export PNETCDFDIR=$HOME/pnetcdf-1.6.1
export HDF5DIR=$HOME/hdf5-1.8.12-intel

export HOMME_ROOT=/ascldap/users/pabosle/ACME/components/homme

rm -rf CMake* Testing/

cmake \
  -D CMAKE_BUILD_TYPE=RELEASE \
  -D CMAKE_C_COMPILER=/opt/openmpi-1.8-intel/bin/mpicc \
  -D CMAKE_CXX_COMPILER=/opt/openmpi-1.8-intel/bin/mpicxx \
  -D CMAKE_Fortran_COMPILER=/opt/openmpi-1.8-intel/bin/mpif90 \
  -D CMAKE_CXX_FLAGS:STRING="-std=c++11" \
  -D HOMME_PROJID:STRING=Aeras \
  -D HOMME_FIND_BLASLAPACK:BOOL=ON \
  -D BUILD_HOMME_SWIM:BOOL=ON \
  -D BUILD_HOMME_PRIM:BOOL=OFF \
  -D BUILD_HOMME_SWEQX:BOOL=ON \
  -D BUILD_HOMME_PREQX:BOOL=OFF \
  -D BUILD_HOMME_SWDGX:BOOL=OFF \
  -D BUILD_HOMME_PRIMDGX:BOOL=OFF \
  -D REFSOLN:BOOL=OFF \
  -D SWIM_PLEV=1 \
  -D SWIM_NP=4 \
  -D SWIM_NC=4 \
  -D PRIM_PLEV=26 \
  -D PRIM_NP=4 \
  -D PRIM_NC=4 \
  -D PIO_FILESYSTEM_HINTS:STRING="lustre" \
  -D CMAKE_INSTALL_PREFIX:PATH=$HOMME_ROOT/build/install \
  -D NETCDF_DIR:PATH=${NETCDFDIR} \
  -D HDF5_DIR:PATH=${HDF5DIR} \
  -D CURL_INCLUDE_DIR:FILEPATH=/usr/bin/curl \
  $HOMME_ROOT