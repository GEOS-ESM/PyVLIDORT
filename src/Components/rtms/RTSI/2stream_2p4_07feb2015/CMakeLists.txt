esma_set_this ()

set (src_dirs .)
#set (srcs LIDORT_Mod.F90 LIDORT_SurfaceMod.F90 LIDORT_ScatMod.F90)
foreach (dir ${src_dirs})
    file (GLOB srcs_ CONFIGURE_DEPENDS ${dir}/*.[Ffc] ${dir}/*.[fF]90)
  list (APPEND srcs ${srcs_})
endforeach()


esma_add_library (${this}
  SRCS ${srcs}
  )
