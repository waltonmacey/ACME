cmake_minimum_required(VERSION 2.7)

OPTION(ACME_BUILD "Flag to distinguish Cmake logic for full ACME build versus Unit Tests" ON)

PROJECT(ACME Fortran C)
FIND_PACKAGE(MPI)

#AGS: needed this or ran into error in  gptl. Need logic?
ADD_DEFINITIONS(-DHAVE_MPI)


# Require installed NETCDF under ACME_NETCDF_DIR
IF (NOT DEFINED ACME_NETCDF_DIR)
  MESSAGE(FATAL_ERROR "ACME_NETCDF_DIR must be defined.")
ENDIF()
IF (EXISTS ${ACME_NETCDF_DIR}/include/netcdf.mod)
  INCLUDE_DIRECTORIES(${ACME_NETCDF_DIR}/include)
  MESSAGE("-- Found netcdf.mod")
ELSE()
  MESSAGE(FATAL_ERROR "Can't find netcdf.mod under ${ACME_NETCDF_DIR}/include")
ENDIF()

# Require installed HDF5 under ACME_HDF5_DIR
IF (NOT DEFINED ACME_HDF5_DIR)
  MESSAGE(FATAL_ERROR "ACME_HDF5_DIR must be defined.")
ENDIF()
IF (EXISTS ${ACME_HDF5_DIR}/lib/libhdf5.a)
  INCLUDE_DIRECTORIES(${ACME_HDF5_DIR}/include)
  MESSAGE("-- Found libhdf5.a")
ELSE()
  MESSAGE(FATAL_ERROR "Can't find libhdf5.a under ${ACME_HDF5_DIR}/include")
ENDIF()

# Require installed PIO under ACME_PIO_DIR
IF (NOT DEFINED ACME_PIO_DIR)
  MESSAGE(FATAL_ERROR "ACME_PIO_DIR must be defined.")
ENDIF()
IF (EXISTS ${ACME_PIO_DIR}/pio.mod)
  INCLUDE_DIRECTORIES(${ACME_PIO_DIR})
  MESSAGE("-- Found pio.mod")
ELSE()
  MESSAGE(FATAL_ERROR "Can't find pio.mod under ${ACME_PIO_DIR}")
ENDIF()

message("\nAGS hack for GNU:  -DCPRGNU") 
ADD_DEFINITIONS(-DCPRGNU)

message("AGS hack2:   -DFORTRANUNDERSCORE") 
ADD_DEFINITIONS(-DFORTRANUNDERSCORE)

# not sure if needed
#message("AGS hack3:   -DHAVE_COMM_F2C -D_MPI") # for MPAS
#ADD_DEFINITIONS(-DHAVE_COMM_F2C -D_MPI)

# Save path to this directory
SET(ACME_MODEL_DIR ${CMAKE_CURRENT_SOURCE_DIR})
SET(CIME_DIR ${ACME_MODEL_DIR}/cime)

SET(MPEU_DIR cime/externals/mct/mpeu)
SET(MCT_DIR cime/externals/mct/mct)
SET(GPTL_DIR cime/share/timing)
SET(ESMF_DIR cime/share/esmf_wrf_timemgr)
SET(SEQ_DIR cime/driver_cpl/shr)
SET(SHR_DIR cime/share/csm_share/shr)
SET(DRIVER_DIR cime/driver_cpl/driver)

SET(MPEU_BINARY_DIR ${ACME_BINARY_DIR}/${MPEU_DIR})
SET(MCT_BINARY_DIR ${ACME_BINARY_DIR}/${MCT_DIR})
SET(GPTL_BINARY_DIR ${ACME_BINARY_DIR}/${GPTL_DIR})
SET(SEQ_BINARY_DIR ${ACME_BINARY_DIR}/${SEQ_DIR})
SET(SHR_BINARY_DIR ${ACME_BINARY_DIR}/${SHR_DIR})
SET(ESMF_BINARY_DIR ${ACME_BINARY_DIR}/${ESMF_DIR})
SET(DRIVER_BINARY_DIR ${ACME_BINARY_DIR}/${DRIVER_DIR})

# Need include from source dir, not build dir like usual
SET(ACME_SHR_DIR ${ACME_MODEL_DIR}/${SHR_DIR}/../include)

MESSAGE("\nConfiguring libraries, utils, share, driver:")

ADD_SUBDIRECTORY(${MPEU_DIR})
ADD_SUBDIRECTORY(${MCT_DIR})
ADD_SUBDIRECTORY(${GPTL_DIR})
ADD_SUBDIRECTORY(${SHR_DIR})
ADD_SUBDIRECTORY(${ESMF_DIR})
ADD_SUBDIRECTORY(${SEQ_DIR})
ADD_SUBDIRECTORY(${DRIVER_DIR})
#
MESSAGE("\nConfiguring Components:")
SET(LND_D_DIR      cime/components/data_comps/dlnd)
SET(LND_S_DIR      cime/components/stub_comps/slnd/cpl)
SET(OCN_D_DIR      cime/components/data_comps/docn)
SET(OCN_S_DIR      cime/components/stub_comps/socn/cpl)
SET(ICE_D_DIR      cime/components/data_comps/dice)
SET(ICE_S_DIR      cime/components/stub_comps/sice/cpl)
SET(GLC_S_DIR      cime/components/stub_comps/sglc/cpl)
SET(ATM_D_DIR      cime/components/data_comps/datm)
SET(ATM_S_DIR      cime/components/stub_comps/satm/cpl)
SET(ROF_D_DIR      cime/components/data_comps/drof)
SET(ROF_S_DIR      cime/components/stub_comps/srof/cpl)
SET(WAV_S_DIR      cime/components/stub_comps/swav/cpl)
SET(LND_D_BINARY_DIR ${ACME_BINARY_DIR}/${LND_D_DIR})
SET(LND_S_BINARY_DIR ${ACME_BINARY_DIR}/${LND_S_DIR})
SET(OCN_D_BINARY_DIR ${ACME_BINARY_DIR}/${OCN_D_DIR})
SET(OCN_S_BINARY_DIR ${ACME_BINARY_DIR}/${OCN_S_DIR})
SET(ICE_D_BINARY_DIR ${ACME_BINARY_DIR}/${ICE_D_DIR})
SET(ICE_S_BINARY_DIR ${ACME_BINARY_DIR}/${ICE_S_DIR})
SET(GLC_S_BINARY_DIR ${ACME_BINARY_DIR}/${GLC_S_DIR})
SET(ATM_D_BINARY_DIR ${ACME_BINARY_DIR}/${ATM_D_DIR})
SET(ATM_S_BINARY_DIR ${ACME_BINARY_DIR}/${ATM_S_DIR})
SET(ROF_D_BINARY_DIR ${ACME_BINARY_DIR}/${ROF_D_DIR})
SET(ROF_S_BINARY_DIR ${ACME_BINARY_DIR}/${ROF_S_DIR})
SET(WAV_S_BINARY_DIR ${ACME_BINARY_DIR}/${WAV_S_DIR})
ADD_SUBDIRECTORY(${LND_D_DIR})
ADD_SUBDIRECTORY(${LND_S_DIR})
ADD_SUBDIRECTORY(${OCN_D_DIR})
ADD_SUBDIRECTORY(${OCN_S_DIR})
ADD_SUBDIRECTORY(${ICE_D_DIR})
ADD_SUBDIRECTORY(${ICE_S_DIR})
ADD_SUBDIRECTORY(${GLC_S_DIR})
ADD_SUBDIRECTORY(${ATM_D_DIR})
ADD_SUBDIRECTORY(${ATM_S_DIR})
ADD_SUBDIRECTORY(${ROF_D_DIR})
ADD_SUBDIRECTORY(${ROF_S_DIR})
ADD_SUBDIRECTORY(${WAV_S_DIR})

## Add MPAS components -- had to run make to generate inc files
##   MPASLI_ONLY code runs but dies at end with mpi error
#SET(MPAS_FRAMEWORK_DIR  components/mpasli/model/src/framework)
#SET(MPAS_EZXML_DIR  components/mpasli/model/src/external/ezxml)
#SET(MPAS_OPERATORS_DIR  components/mpasli/model/src/operators)
#SET(MPAS_GLC_DIR  components/mpasli/model/src/core_landice)
#SET(MPAS_GLCDRIVER_DIR  components/mpasli/driver)
#SET(MPAS_FRAMEWORK_BINARY_DIR ${ACME_BINARY_DIR}/${MPAS_FRAMEWORK_DIR})
#SET(MPAS_EZXML_BINARY_DIR ${ACME_BINARY_DIR}/${MPAS_EZXML_DIR})
#SET(MPAS_OPERATORS_BINARY_DIR ${ACME_BINARY_DIR}/${MPAS_OPERATORS_DIR})
#SET(MPAS_GLC_BINARY_DIR ${ACME_BINARY_DIR}/${MPAS_GLC_DIR})
#SET(MPAS_GLCDRIVER_BINARY_DIR ${ACME_BINARY_DIR}/${MPAS_GLCDRIVER_DIR})
#
#SET(ACME_MPAS_EZXML_DIR ${ACME_MODEL_DIR}/${MPAS_EZXML_DIR}) # path to ezxml.h
#
#ADD_SUBDIRECTORY(${MPAS_FRAMEWORK_DIR})
#ADD_SUBDIRECTORY(${MPAS_EZXML_DIR})
#ADD_SUBDIRECTORY(${MPAS_OPERATORS_DIR})
#ADD_SUBDIRECTORY(${MPAS_GLC_DIR})
#ADD_SUBDIRECTORY(${MPAS_GLCDRIVER_DIR})

SET(LND_CLM_DIR      components/clm/src)
SET(LND_CLM_BINARY_DIR ${ACME_BINARY_DIR}/${LND_CLM_DIR})
ADD_SUBDIRECTORY(${LND_CLM_DIR})

SET(RTM_DIR      components/rtm/src)
SET(RTM_BINARY_DIR ${ACME_BINARY_DIR}/${RTM_DIR})
ADD_SUBDIRECTORY(${RTM_DIR})


MESSAGE("\nAdding executable:")
MESSAGE("   + acme_A.exe: A-Compset")
ADD_EXECUTABLE(acme_A.exe ${DRIVER_DIR}/cesm_driver.F90)
#ADD_EXECUTABLE(acme_A.exe ${DRIVER_DIR}/cesm_driver.F90 ${DRIVER_DIR}/cesm_comp_mod.F90)
MESSAGE("   + acme_S.exe: S-Compset")
ADD_EXECUTABLE(acme_S.exe ${DRIVER_DIR}/cesm_driver.F90)
#ADD_EXECUTABLE(acme_S.exe ${DRIVER_DIR}/cesm_driver.F90 ${DRIVER_DIR}/cesm_comp_mod.F90)
MESSAGE("   + acme_I.exe: I-Compset")
ADD_EXECUTABLE(acme_I.exe ${DRIVER_DIR}/cesm_driver.F90)
#ADD_EXECUTABLE(acme_I.exe ${DRIVER_DIR}/cesm_driver.F90 ${DRIVER_DIR}/cesm_comp_mod.F90)
#MESSAGE("   + acme_MPASLI_ONLY.exe: MPASLI_ONLY-Compset")
#ADD_EXECUTABLE(acme_MPASLI_ONLY.exe ${DRIVER_DIR}/cesm_driver.F90 ${DRIVER_DIR}/cesm_comp_mod.F90)


MESSAGE("\nAGS: HHEELLPP  what include dirs to add: data, stub, active?
Executable dependent??")

INCLUDE_DIRECTORIES (${SHR_BINARY_DIR}
                     ${GPTL_BINARY_DIR}
                     ${ESMF_BINARY_DIR}
                     ${SEQ_BINARY_DIR}
                     ${DRIVER_BINARY_DIR}
                     ${LND_D_BINARY_DIR}
                     ${OCN_D_BINARY_DIR}
                     ${ICE_D_BINARY_DIR}
                     ${GLC_S_BINARY_DIR}
                     ${ATM_D_BINARY_DIR}
                     ${ROF_D_BINARY_DIR}
                     ${WAV_S_BINARY_DIR})

MESSAGE("   + CompMod:")
ADD_LIBRARY(compmod ${DRIVER_DIR}/cesm_comp_mod.F90)
ADD_DEPENDENCIES(compmod dlnd docn dice sglc datm drof swav driver)

LINK_DIRECTORIES( ${ACME_NETCDF_DIR}/lib ${ACME_PIO_DIR})
set (PIO_LIB ${ACME_PIO_DIR}/libpio.a)
set (NETCDF_LIB ${ACME_NETCDF_DIR}/lib/libnetcdff.a
                ${ACME_NETCDF_DIR}/lib/libnetcdf.a
                ${ACME_NETCDF_DIR}/lib/libpnetcdf.a
                ${ACME_HDF5_DIR}/lib/libhdf5_hl.a
                ${ACME_HDF5_DIR}/lib/libhdf5.a
                z )

message ("\nAGS: these are pio netcdf Libs : ${PIO_LIB} ${NETCDF_LIB}")

TARGET_LINK_LIBRARIES(acme_A.exe  compmod
                      dlnd docn dice sglc datm drof swav 
                      driver shr seq shrcal gptl esmf mct mpeu shrsys 
                      ${PIO_LIB} ${NETCDF_LIB})

TARGET_LINK_LIBRARIES(acme_S.exe  compmod
                      slnd socn sice sglc satm srof swav 
                      driver shr seq shrcal gptl esmf mct mpeu shrsys 
                      ${PIO_LIB} ${NETCDF_LIB})

TARGET_LINK_LIBRARIES(acme_I.exe  compmod
                      clm socn sice sglc datm rtm swav 
                      driver shr seq shrcal gptl esmf mct mpeu shrsys 
                      ${PIO_LIB} ${NETCDF_LIB} -llapack)

#TARGET_LINK_LIBRARIES(acme_MPASLI_ONLY.exe  
#                      slnd socn sice 
#                        mpasglcdriver mpasglc mpasoperators mpasframework mpasezxml
#                      satm srof swav 
#                      driver shr seq shrcal gptl esmf mct mpeu shrsys 
#                      ${PIO_LIB} ${NETCDF_LIB})

MESSAGE("\nAnd CMake says...:")
