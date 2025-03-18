# We only want to build and use one VLIDORT library

set(VLIDORT_DEFAULT 2p8p3)
set(VLIDORT ${VLIDORT_DEFAULT} CACHE STRING "VLIDORT Library to link to")
string(CONCAT VLIDORT_LIBRARY "VLIDORT" ${VLIDORT})


f(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${VLIDORT_LIBRARY}" AND IS_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/${VLIDORT_LIBRARY}")
  message(STATUS "Building with VLIDORT library ${VLIDORT_LIBRARY}")
else()
  message(FATAL_ERROR "VLIDORT library ${VLIDORT_LIBRARY} does not exist")
endif()

esma_add_subdirectories (2STREAM LIDORT3p7 ${VLIDORT_LIBRARY})
