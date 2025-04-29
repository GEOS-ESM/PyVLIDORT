esma_set_this ()

set (src_dirs 2stream_sourcecode)
set (srcs TWOSTREAM_Mod.F90 TWOSTREAM_SurfaceMod.F90 TWOSTREAM_ScatMod.F90 TWOSTREAM_pars.F90)
foreach (dir ${src_dirs})
    file (GLOB srcs_ CONFIGURE_DEPENDS ${dir}/*.[Ffc] ${dir}/*.[fF]90)
  list (APPEND srcs ${srcs_})
endforeach()


esma_add_library (${this}
  SRCS ${srcs}
  )
