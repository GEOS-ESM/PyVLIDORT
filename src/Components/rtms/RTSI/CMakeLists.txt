
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${VLIDORT_LIBRARY}" AND IS_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/${VLIDORT_LIBRARY}")
  message(STATUS "Building with VLIDORT library ${VLIDORT_LIBRARY}")
else()
  message(FATAL_ERROR "VLIDORT library ${VLIDORT_LIBRARY} does not exist")
endif()

esma_add_subdirectories (LIDORT3p7 ${VLIDORT_LIBRARY} 2stream_2p4_TDCSModified_08Apr24)
