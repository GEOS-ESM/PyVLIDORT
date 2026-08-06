esma_set_this ()

set (src_dirs 2stream_sourcecode)
set (srcs TWOSTREAM_Mod.F90 TWOSTREAM_SurfaceMod.F90 TWOSTREAM_ScatMod.F90 TWOSTREAM_pars.F90 
            twostream_master_debug_Mod.F90 twostream_master_outputs_debug_Mod.F90
            vfzmat_Master_debug_mod.f90)
foreach (dir ${src_dirs})
    file (GLOB srcs_ CONFIGURE_DEPENDS ${dir}/*.[Ffc] ${dir}/*.[fF]90)
  list (APPEND srcs ${srcs_})
endforeach()

set(VFZMAT_PATH vfzmat)
list (APPEND srcs
   ${VFZMAT_PATH}/vfzmat_Numerical.f90
   ${VFZMAT_PATH}/vfzmat_DevelopCoeffs.f90
   ${VFZMAT_PATH}/vfzmat_ExpandCoeffs.f90
   ${VFZMAT_PATH}/vfzmat_Rotation.f90
   ${VFZMAT_PATH}/vfzmat_PhasMat.f90
   ${VFZMAT_PATH}/vfzmat_Rayleigh.f90
   ${VFZMAT_PATH}/vfzmat_Master.f90
)

esma_add_library (${this}
  SRCS ${srcs}
  )
