MESSAGE("   + MPAS glcdriver:")
FILE(GLOB MPAS_GLCDRIVER_SOURCE *F)

INCLUDE_DIRECTORIES (
      ${SEQ_BINARY_DIR}
      ${SHR_BINARY_DIR}
      ${ESMF_BINARY_DIR}
      ${GPTL_BINARY_DIR}
      ${MPAS_FRAMEWORK_BINARY_DIR}
      ${MPAS_GLC_BINARY_DIR}
)

ADD_DEFINITIONS(-ffree-form)
 
ADD_LIBRARY (mpasglcdriver ${MPAS_GLCDRIVER_SOURCE})
ADD_DEPENDENCIES(mpasglcdriver  seq shrcal esmf gptl mpasframework mpasglc)
