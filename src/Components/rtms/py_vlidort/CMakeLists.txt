esma_set_this ()
find_package(F2PY3 REQUIRED)


# MAT: Note uses add_f2py3 because I can't figure out the import test for netCDF f2py code
add_f2py3_module (VLIDORT_POLAR_
  SOURCES VLIDORT_POLAR_py.F90
  DESTINATION lib/Python/${this}
  LIBRARIES drivers ${VLIDORT_LIBRARY}
  INCLUDEDIRS ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_BINARY_DIR}/lib ${include_${this}} ${esma_include}/drivers ${esma_include}/${VLIDORT_LIBRARY}
  )
add_dependencies(VLIDORT_POLAR_ drivers ${VLIDORT_LIBRARY})

add_f2py3_module (TWOSTREAM_
  SOURCES TWOSTREAM_py.F90
  DESTINATION lib/Python/${this}
  LIBRARIES drivers TWOSTREAM2p4
  INCLUDEDIRS ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_BINARY_DIR}/lib ${include_${this}} ${esma_include}/drivers ${esma_include}/TWOSTREAM2p4
  )
add_dependencies(TWOSTREAM_ drivers TWOSTREAM2p4)

file (GLOB python_files CONFIGURE_DEPENDS py_vlidort/*.py)
install (FILES ${python_files} DESTINATION lib/Python/py_vlidort)

