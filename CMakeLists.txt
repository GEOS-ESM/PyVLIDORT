#
# Top cmake rules for pyVLIDORT
#

# Boiler plate preamble
# ---------------------
cmake_minimum_required (VERSION 3.13)
cmake_policy (SET CMP0053 NEW)
cmake_policy (SET CMP0054 NEW)



esma_set_this()

set (alldirs
  )

esma_add_library (${this}
  SRCS 
  SUBCOMPONENTS ${alldirs}
  DEPENDENCIES MAPL esmf)

