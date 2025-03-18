# We only want to build and use one VLIDORT library

set(VLIDORT_DEFAULT 2p8p3)
set(VLIDORT ${VLIDORT_DEFAULT} CACHE STRING "VLIDORT Library to link to")
string(CONCAT VLIDORT_LIBRARY "VLIDORT" ${VLIDORT})


add_subdirectory (drivers)
add_subdirectory (py_vlidort)
add_subdirectory (RTSI)
