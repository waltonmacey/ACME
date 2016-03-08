# CMake initial cache file for Cori Phase 1 (corip1)

SET (CMAKE_Fortran_COMPILER ftn CACHE FILEPATH "")
SET (CMAKE_C_COMPILER cc CACHE FILEPATH "")
SET (CMAKE_CXX_COMPILER CC CACHE FILEPATH "")

SET (NETCDF_DIR $ENV{NETCDF_DIR} CACHE FILEPATH "")
#example with module cray-netcdf-hdf5parallel/4.3.3.1:
# NETCDF_DIR=/opt/cray/netcdf-hdf5parallel/4.3.3.1/INTEL/14.0

#ndk SET (PNETCDF_DIR $ENV{PARALLEL_NETCDF_DIR} CACHE FILEPATH "")
# this env var is not set with module cray-netcdf-hdf5parallel/4.3.3.1

SET (HDF5_DIR $ENV{HDF5_DIR} CACHE FILEPATH "")
#example with module cray-hdf5-parallel/1.8.14: 
# HDF5_DIR=/opt/cray/hdf5-parallel/1.8.14/INTEL/14.0

SET (CMAKE_SYSTEM_NAME Catamount CACHE FILEPATH "")
