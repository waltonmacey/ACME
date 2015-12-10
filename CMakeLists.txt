cmake_minimum_required(VERSION 2.7)

OPTION(ACME_BUILD "Flag to distinguish Cmake logic for full ACME build versus Unit Tests" ON)

PROJECT(ACME Fortran C)
FIND_PACKAGE(MPI)


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
SET(LND_DIR      cime/components/data_comps/dlnd)
SET(OCN_DIR      cime/components/data_comps/docn)
SET(ICE_DIR      cime/components/data_comps/dice)
SET(GLC_DIR      cime/components/stub_comps/sglc/cpl)
SET(ATM_DIR      cime/components/data_comps/datm)
SET(ROF_DIR      cime/components/data_comps/drof)
SET(WAV_DIR      cime/components/stub_comps/swav/cpl)
SET(LND_BINARY_DIR ${ACME_BINARY_DIR}/${LND_DIR})
SET(OCN_BINARY_DIR ${ACME_BINARY_DIR}/${OCN_DIR})
SET(ICE_BINARY_DIR ${ACME_BINARY_DIR}/${ICE_DIR})
SET(GLC_BINARY_DIR ${ACME_BINARY_DIR}/${GLC_DIR})
SET(ATM_BINARY_DIR ${ACME_BINARY_DIR}/${ATM_DIR})
SET(ROF_BINARY_DIR ${ACME_BINARY_DIR}/${ROF_DIR})
SET(WAV_BINARY_DIR ${ACME_BINARY_DIR}/${WAV_DIR})
ADD_SUBDIRECTORY(${LND_DIR})
ADD_SUBDIRECTORY(${OCN_DIR})
ADD_SUBDIRECTORY(${ICE_DIR})
ADD_SUBDIRECTORY(${GLC_DIR})
ADD_SUBDIRECTORY(${ATM_DIR})
ADD_SUBDIRECTORY(${ROF_DIR})
ADD_SUBDIRECTORY(${WAV_DIR})

# ToDo: add active components
#ADD_SUBDIRECTORY(components/clm/src)
#ADD_SUBDIRECTORY(components/cam)

MESSAGE("\nAdding executable:")
MESSAGE("   + acme.exe:")
# Need to compile drv/driver/ccsm_driver.F90 ccsm_comp_mod.F90, and link
ADD_EXECUTABLE(acme.exe
    ${DRIVER_DIR}/cesm_driver.F90 ${DRIVER_DIR}/cesm_comp_mod.F90
)
INCLUDE_DIRECTORIES (${SHR_BINARY_DIR}
                     ${GPTL_BINARY_DIR}
                     ${ESMF_BINARY_DIR}
                     ${SEQ_BINARY_DIR}
                     ${DRIVER_BINARY_DIR}
                     ${LND_BINARY_DIR}
                     ${OCN_BINARY_DIR}
                     ${ICE_BINARY_DIR}
                     ${GLC_BINARY_DIR}
                     ${ATM_BINARY_DIR}
                     ${ROF_BINARY_DIR}
                     ${WAV_BINARY_DIR})

LINK_DIRECTORIES( ${ACME_NETCDF_DIR}/lib ${ACME_PIO_DIR})
set (PIO_LIB ${ACME_PIO_DIR}/libpio.a)
set (NETCDF_LIB ${ACME_NETCDF_DIR}/lib/libnetcdff.a
                ${ACME_NETCDF_DIR}/lib/libnetcdf.a
                ${ACME_NETCDF_DIR}/lib/libpnetcdf.a
                /home/agsalin/install/hdf5-1.8.11/hdf5/lib/libhdf5_hl.a
                /home/agsalin/install/hdf5-1.8.11/hdf5/lib/libhdf5.a
                z )

message ("\nAGS: these are pio netcdf Libs : ${PIO_LIB} ${NETCDF_LIB}")

TARGET_LINK_LIBRARIES(acme.exe  
                      dlnd docn dice sglc datm drof swav 
                      driver shr seq shrcal gptl esmf mct mpeu shrsys 
                      ${PIO_LIB} ${NETCDF_LIB})

MESSAGE("\nAnd CMake says...:")
