module vfzmat_Rayleigh_debug_m

  implicit none

contains


subroutine write_vfzmat_Rayleigh_debug( &
    MAX_GEOMETRIES, MAX_SZANGLES, MAX_USER_VZANGLES, MAX_USER_RELAZMS, &
    MAXLAYERS,                                                          &
    do_upwelling, do_dnwelling, do_ObsGeoms, do_Sunlight,              &
    NSTOKES, N_GEOMS, N_SZAS, N_VZAS, N_AZMS,                         &
    OFFSETS, DEG_TO_RAD, SZAS, VZAS, AZMS, OBSGEOMS, DEPOL_RATIO,     &
    RayFmatrices_up, RayFmatrices_dn,                                  &
    RayZmatrices_up, RayZmatrices_dn,                                  &
    RayCoeffs )

  implicit none

  !--- Input arguments ---
  integer, intent(in) :: MAX_GEOMETRIES, MAX_SZANGLES
  integer, intent(in) :: MAX_USER_VZANGLES, MAX_USER_RELAZMS
  integer, intent(in) :: MAXLAYERS
  logical, intent(in) :: do_upwelling, do_dnwelling
  logical, intent(in) :: do_ObsGeoms, do_Sunlight
  integer, intent(in) :: NSTOKES, N_GEOMS, N_SZAS, N_VZAS, N_AZMS
  integer, intent(in) :: OFFSETS(MAX_SZANGLES, MAX_USER_VZANGLES)
  real*8,  intent(in) :: DEG_TO_RAD
  real*8,  intent(in) :: SZAS(MAX_SZANGLES)
  real*8,  intent(in) :: VZAS(MAX_USER_VZANGLES)
  real*8,  intent(in) :: AZMS(MAX_USER_RELAZMS)
  real*8,  intent(in) :: OBSGEOMS(MAX_GEOMETRIES, 3)
  real*8,  intent(in) :: DEPOL_RATIO
  real*8,  intent(in) :: RayFmatrices_up(MAX_GEOMETRIES, 6)
  real*8,  intent(in) :: RayFmatrices_dn(MAX_GEOMETRIES, 6)
  real*8,  intent(in) :: RayZmatrices_up(MAX_GEOMETRIES, 4, 4)
  real*8,  intent(in) :: RayZmatrices_dn(MAX_GEOMETRIES, 4, 4)
  real*8,  intent(in) :: RayCoeffs(0:2, 6)

  !--- Local variables ---
  integer :: i, j, k

  !--- Open output file ---
  OPEN(UNIT=99, FILE='vfzmat_Rayleigh_debug.txt', STATUS='REPLACE', ACTION='WRITE')

  !========================================================================
  WRITE(99,'(A)') '========================================================'
  WRITE(99,'(A)') '   vfzmat_Rayleigh: INPUT/OUTPUT DIAGNOSTIC DUMP'
  WRITE(99,'(A)') '========================================================'

  !------------------------------------------------------------------------
  WRITE(99,'(A)') ''
  WRITE(99,'(A)') '--- INPUT DIMENSIONS ---'
  WRITE(99,'(A,I0)') '  MAX_GEOMETRIES       = ', MAX_GEOMETRIES
  WRITE(99,'(A,I0)') '  MAX_SZANGLES         = ', MAX_SZANGLES
  WRITE(99,'(A,I0)') '  MAX_USER_VZANGLES    = ', MAX_USER_VZANGLES
  WRITE(99,'(A,I0)') '  MAX_USER_RELAZMS     = ', MAX_USER_RELAZMS
  WRITE(99,'(A,I0)') '  MAXLAYERS            = ', MAXLAYERS

  !------------------------------------------------------------------------
  WRITE(99,'(A)') ''
  WRITE(99,'(A)') '--- INPUT FLAGS ---'
  WRITE(99,'(A,L1)') '  do_upwelling         = ', do_upwelling
  WRITE(99,'(A,L1)') '  do_dnwelling         = ', do_dnwelling
  WRITE(99,'(A,L1)') '  do_ObsGeoms          = ', do_ObsGeoms
  WRITE(99,'(A,L1)') '  do_Sunlight          = ', do_Sunlight

  !------------------------------------------------------------------------
  WRITE(99,'(A)') ''
  WRITE(99,'(A)') '--- INPUT NUMBERS ---'
  WRITE(99,'(A,I0)') '  NSTOKES              = ', NSTOKES
  WRITE(99,'(A,I0)') '  N_GEOMS              = ', N_GEOMS
  WRITE(99,'(A,I0)') '  N_SZAS               = ', N_SZAS
  WRITE(99,'(A,I0)') '  N_VZAS               = ', N_VZAS
  WRITE(99,'(A,I0)') '  N_AZMS               = ', N_AZMS

  !------------------------------------------------------------------------
  WRITE(99,'(A)') ''
  WRITE(99,'(A)') '--- INPUT SCALARS ---'
  WRITE(99,'(A,E20.10)') '  DEG_TO_RAD           = ', DEG_TO_RAD
  WRITE(99,'(A,E20.10)') '  DEPOL_RATIO          = ', DEPOL_RATIO

  !------------------------------------------------------------------------
  WRITE(99,'(A)') ''
  WRITE(99,'(A)') '--- INPUT OFFSETS (N_SZAS x N_VZAS active values) ---'
  WRITE(99,'(A)') '  Row = SZA index, Col = VZA index'
  DO i = 1, N_SZAS
    WRITE(99,'(A,I4,A)', ADVANCE='NO') '  SZA(', i, '): '
    DO j = 1, N_VZAS
      WRITE(99,'(I6)', ADVANCE='NO') OFFSETS(i,j)
    END DO
    WRITE(99,*)
  END DO

  !------------------------------------------------------------------------
  WRITE(99,'(A)') ''
  WRITE(99,'(A)') '--- INPUT SZAS (degrees) ---'
  DO i = 1, N_SZAS
    WRITE(99,'(A,I4,A,F12.6)') '  SZAS(', i, ') = ', SZAS(i)
  END DO

  !------------------------------------------------------------------------
  WRITE(99,'(A)') ''
  WRITE(99,'(A)') '--- INPUT VZAS (degrees) ---'
  DO i = 1, N_VZAS
    WRITE(99,'(A,I4,A,F12.6)') '  VZAS(', i, ') = ', VZAS(i)
  END DO

  !------------------------------------------------------------------------
  WRITE(99,'(A)') ''
  WRITE(99,'(A)') '--- INPUT AZMS (degrees) ---'
  DO i = 1, N_AZMS
    WRITE(99,'(A,I4,A,F12.6)') '  AZMS(', i, ') = ', AZMS(i)
  END DO

  !------------------------------------------------------------------------
  WRITE(99,'(A)') ''
  WRITE(99,'(A)') '--- INPUT OBSGEOMS (N_GEOMS active rows x 3) ---'
  WRITE(99,'(A)') '  Columns: [SZA, VZA, AZM]'
  DO i = 1, N_GEOMS
    WRITE(99,'(A,I4,A,3F12.6)') '  OBSGEOMS(', i, ',:) = ', &
      OBSGEOMS(i,1), OBSGEOMS(i,2), OBSGEOMS(i,3)
  END DO

  !========================================================================
  WRITE(99,'(A)') ''
  WRITE(99,'(A)') '========================================================'
  WRITE(99,'(A)') '--- OUTPUTS ---'
  WRITE(99,'(A)') '========================================================'

  !------------------------------------------------------------------------
  WRITE(99,'(A)') ''
  WRITE(99,'(A)') '--- OUTPUT RayFmatrices_up (N_GEOMS x 6) ---'
  DO i = 1, N_GEOMS
      WRITE(99,'(A,I4,A,6E16.6)') &
        '  geom=', i, '  : ', RayFmatrices_up(i,1:6)
  END DO

  !------------------------------------------------------------------------
  WRITE(99,'(A)') ''
  WRITE(99,'(A)') '--- OUTPUT RayFmatrices_dn (N_GEOMS x MAXLAYERS x 6) ---'
  DO i = 1, N_GEOMS
      WRITE(99,'(A,I4,A,6E16.6)') &
        '  geom=', i, '  : ', RayFmatrices_dn(i,1:6)
  END DO

  !------------------------------------------------------------------------
  WRITE(99,'(A)') ''
  WRITE(99,'(A)') '--- OUTPUT RayZmatrices_up (N_GEOMS x 4 x 4) ---'
  DO i = 1, N_GEOMS
    WRITE(99,'(A,I4)') '  geom=', i
    DO k = 1, 4
      WRITE(99,'(A,4E16.6)') '    row: ', RayZmatrices_up(i,k,1:4)
    END DO
  END DO

  !------------------------------------------------------------------------
  WRITE(99,'(A)') ''
  WRITE(99,'(A)') '--- OUTPUT RayZmatrices_dn (N_GEOMS x 4 x 4) ---'
  DO i = 1, N_GEOMS
    WRITE(99,'(A,I4)') '  geom=', i
    DO k = 1, 4
      WRITE(99,'(A,4E16.6)') '    row: ', RayZmatrices_dn(i,k,1:4)
    END DO
  END DO

  !------------------------------------------------------------------------
  WRITE(99,'(A)') ''
  WRITE(99,'(A)') '--- OUTPUT RayCoeffs (0:2 x 6) ---'
  DO k = 0, 2
    WRITE(99,'(A,I2,A,6E16.6)') '  moment=', k, '  : ', RayCoeffs(k,1:6)
  END DO

  !========================================================================
  WRITE(99,'(A)') ''
  WRITE(99,'(A)') '========================================================'
  WRITE(99,'(A)') '  END OF DUMP'
  WRITE(99,'(A)') '========================================================'

  CLOSE(99)

end subroutine write_vfzmat_Rayleigh_debug

end module vfzmat_Rayleigh_debug_m
