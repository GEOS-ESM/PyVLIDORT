module vfzmat_Master_debug_m

  implicit none

contains

  subroutine write_vfzmat_Master_debug( &
      MAXMOMENTS_INPUT, MAX_GEOMETRIES, MAX_SZANGLES, MAX_USER_VZANGLES, &
      MAX_USER_RELAZMS, MAXLAYERS,                                        &
      Max_InAngles, N_InAngles, InAngles, fmatrix, Exist_InFmatrices,    &
      do_upwelling, do_dnwelling, do_ObsGeoms, do_Sunlight,              &
      nmom, NLAYERS, NSTOKES, N_GEOMS, N_SZAS, N_VZAS, N_AZMS,          &
      OFFSETS, DEG_TO_RAD, SZAS, VZAS, AZMS, OBSGEOMS,                  &
      AerFmatrices_up, AerFmatrices_dn,                                  &
      AerZmatrices_up, AerZmatrices_dn,                                  &
      AerCoeffs )

    implicit none

    !--- Input arguments ---
    integer, intent(in) :: MAXMOMENTS_INPUT, MAX_GEOMETRIES, MAX_SZANGLES
    integer, intent(in) :: MAX_USER_VZANGLES, MAX_USER_RELAZMS, MAXLAYERS
    integer, intent(in) :: Max_InAngles, N_InAngles
    real*8,  intent(in) :: InAngles(Max_InAngles)
    real*8,  intent(in) :: fmatrix(NLAYERS, 1, 1, N_InAngles, 6)
    logical, intent(in) :: Exist_InFmatrices(MAXLAYERS)
    logical, intent(in) :: do_upwelling, do_dnwelling
    logical, intent(in) :: do_ObsGeoms, do_Sunlight
    integer, intent(in) :: nmom, NLAYERS, NSTOKES
    integer, intent(in) :: N_GEOMS, N_SZAS, N_VZAS, N_AZMS
    integer, intent(in) :: OFFSETS(MAX_SZANGLES, MAX_USER_VZANGLES)
    real*8,  intent(in) :: DEG_TO_RAD
    real*8,  intent(in) :: SZAS(MAX_SZANGLES)
    real*8,  intent(in) :: VZAS(MAX_USER_VZANGLES)
    real*8,  intent(in) :: AZMS(MAX_USER_RELAZMS)
    real*8,  intent(in) :: OBSGEOMS(MAX_GEOMETRIES, 3)
    real*8,  intent(in) :: AerFmatrices_up(MAX_GEOMETRIES, MAXLAYERS, 6)
    real*8,  intent(in) :: AerFmatrices_dn(MAX_GEOMETRIES, MAXLAYERS, 6)
    real*8,  intent(in) :: AerZmatrices_up(MAX_GEOMETRIES, MAXLAYERS, 4, 4)
    real*8,  intent(in) :: AerZmatrices_dn(MAX_GEOMETRIES, MAXLAYERS, 4, 4)
    real*8,  intent(in) :: AerCoeffs(MAXLAYERS, 0:MAXMOMENTS_INPUT, 6)

    !--- Local variables ---
    integer :: i, j, k, m, p

    !--- Open output file ---
    OPEN(UNIT=98, FILE='vfzmat_Master_debug.txt', STATUS='REPLACE', ACTION='WRITE')

    !========================================================================
    WRITE(98,'(A)') '========================================================'
    WRITE(98,'(A)') '   vfzmat_Master: INPUT/OUTPUT DIAGNOSTIC DUMP'
    WRITE(98,'(A)') '========================================================'

    !------------------------------------------------------------------------
    WRITE(98,'(A)') ''
    WRITE(98,'(A)') '--- INPUT DIMENSIONS ---'
    WRITE(98,'(A,I0)') '  MAXMOMENTS_INPUT     = ', MAXMOMENTS_INPUT
    WRITE(98,'(A,I0)') '  MAX_GEOMETRIES       = ', MAX_GEOMETRIES
    WRITE(98,'(A,I0)') '  MAX_SZANGLES         = ', MAX_SZANGLES
    WRITE(98,'(A,I0)') '  MAX_USER_VZANGLES    = ', MAX_USER_VZANGLES
    WRITE(98,'(A,I0)') '  MAX_USER_RELAZMS     = ', MAX_USER_RELAZMS
    WRITE(98,'(A,I0)') '  MAXLAYERS            = ', MAXLAYERS

    !------------------------------------------------------------------------
    WRITE(98,'(A)') ''
    WRITE(98,'(A)') '--- INPUT FMATRIX INFO ---'
    WRITE(98,'(A,I0)') '  Max_InAngles         = ', Max_InAngles
    WRITE(98,'(A,I0)') '  N_InAngles           = ', N_InAngles

    WRITE(98,'(A)') ''
    WRITE(98,'(A)') '  InAngles (active values):'
    DO i = 1, N_InAngles
      WRITE(98,'(A,I4,A,F12.6)') '    InAngles(', i, ') = ', InAngles(i)
    END DO

    WRITE(98,'(A)') ''
    WRITE(98,'(A)') '  Exist_InFmatrices (layer flags):'
    DO k = 1, NLAYERS
      WRITE(98,'(A,I4,A,L1)') '    layer=', k, '  : ', Exist_InFmatrices(k)
    END DO

    WRITE(98,'(A)') ''
    WRITE(98,'(A)') '  fmatrix (NLAYERS x nch=1 x nobs=1 x N_InAngles x nPol=6):'
    WRITE(98,'(A)') '  [layer, angle, pol]'
    DO k = 1, NLAYERS
      DO i = 1, N_InAngles
        WRITE(98,'(A,I4,A,I4,A,6E16.6)') &
          '    layer=', k, '  angle=', i, '  : ', fmatrix(k,1,1,i,1:6)
      END DO
    END DO

    !------------------------------------------------------------------------
    WRITE(98,'(A)') ''
    WRITE(98,'(A)') '--- INPUT FLAGS ---'
    WRITE(98,'(A,L1)') '  do_upwelling         = ', do_upwelling
    WRITE(98,'(A,L1)') '  do_dnwelling         = ', do_dnwelling
    WRITE(98,'(A,L1)') '  do_ObsGeoms          = ', do_ObsGeoms
    WRITE(98,'(A,L1)') '  do_Sunlight          = ', do_Sunlight

    !------------------------------------------------------------------------
    WRITE(98,'(A)') ''
    WRITE(98,'(A)') '--- INPUT NUMBERS ---'
    WRITE(98,'(A,I0)') '  nmom                 = ', nmom
    WRITE(98,'(A,I0)') '  NLAYERS              = ', NLAYERS
    WRITE(98,'(A,I0)') '  NSTOKES              = ', NSTOKES
    WRITE(98,'(A,I0)') '  N_GEOMS              = ', N_GEOMS
    WRITE(98,'(A,I0)') '  N_SZAS               = ', N_SZAS
    WRITE(98,'(A,I0)') '  N_VZAS               = ', N_VZAS
    WRITE(98,'(A,I0)') '  N_AZMS               = ', N_AZMS

    !------------------------------------------------------------------------
    WRITE(98,'(A)') ''
    WRITE(98,'(A)') '--- INPUT SCALARS ---'
    WRITE(98,'(A,E20.10)') '  DEG_TO_RAD           = ', DEG_TO_RAD

    !------------------------------------------------------------------------
    WRITE(98,'(A)') ''
    WRITE(98,'(A)') '--- INPUT OFFSETS (N_SZAS x N_VZAS active values) ---'
    WRITE(98,'(A)') '  Row = SZA index, Col = VZA index'
    DO i = 1, N_SZAS
      WRITE(98,'(A,I4,A)', ADVANCE='NO') '  SZA(', i, '): '
      DO j = 1, N_VZAS
        WRITE(98,'(I6)', ADVANCE='NO') OFFSETS(i,j)
      END DO
      WRITE(98,*)
    END DO

    !------------------------------------------------------------------------
    WRITE(98,'(A)') ''
    WRITE(98,'(A)') '--- INPUT SZAS (degrees) ---'
    DO i = 1, N_SZAS
      WRITE(98,'(A,I4,A,F12.6)') '  SZAS(', i, ') = ', SZAS(i)
    END DO

    !------------------------------------------------------------------------
    WRITE(98,'(A)') ''
    WRITE(98,'(A)') '--- INPUT VZAS (degrees) ---'
    DO i = 1, N_VZAS
      WRITE(98,'(A,I4,A,F12.6)') '  VZAS(', i, ') = ', VZAS(i)
    END DO

    !------------------------------------------------------------------------
    WRITE(98,'(A)') ''
    WRITE(98,'(A)') '--- INPUT AZMS (degrees) ---'
    DO i = 1, N_AZMS
      WRITE(98,'(A,I4,A,F12.6)') '  AZMS(', i, ') = ', AZMS(i)
    END DO

    !------------------------------------------------------------------------
    WRITE(98,'(A)') ''
    WRITE(98,'(A)') '--- INPUT OBSGEOMS (N_GEOMS active rows x 3) ---'
    WRITE(98,'(A)') '  Columns: [SZA, VZA, AZM]'
    DO i = 1, N_GEOMS
      WRITE(98,'(A,I4,A,3F12.6)') '  OBSGEOMS(', i, ',:) = ', &
        OBSGEOMS(i,1), OBSGEOMS(i,2), OBSGEOMS(i,3)
    END DO

    !========================================================================
    WRITE(98,'(A)') ''
    WRITE(98,'(A)') '========================================================'
    WRITE(98,'(A)') '--- OUTPUTS ---'
    WRITE(98,'(A)') '========================================================'

    !------------------------------------------------------------------------
    WRITE(98,'(A)') ''
    WRITE(98,'(A)') '--- OUTPUT AerFmatrices_up (N_GEOMS x NLAYERS x 6) ---'
    DO i = 1, N_GEOMS
      DO k = 1, NLAYERS
        WRITE(98,'(A,I4,A,I4,A,6E16.6)') &
          '  geom=', i, '  layer=', k, '  : ', AerFmatrices_up(i,k,1:6)
      END DO
    END DO

    !------------------------------------------------------------------------
    WRITE(98,'(A)') ''
    WRITE(98,'(A)') '--- OUTPUT AerFmatrices_dn (N_GEOMS x NLAYERS x 6) ---'
    DO i = 1, N_GEOMS
      DO k = 1, NLAYERS
        WRITE(98,'(A,I4,A,I4,A,6E16.6)') &
          '  geom=', i, '  layer=', k, '  : ', AerFmatrices_dn(i,k,1:6)
      END DO
    END DO

    !------------------------------------------------------------------------
    WRITE(98,'(A)') ''
    WRITE(98,'(A)') '--- OUTPUT AerZmatrices_up (N_GEOMS x NLAYERS x 4 x 4) ---'
    DO i = 1, N_GEOMS
      DO k = 1, NLAYERS
        WRITE(98,'(A,I4,A,I4)') '  geom=', i, '  layer=', k
        DO m = 1, 4
          WRITE(98,'(A,4E16.6)') '    row: ', AerZmatrices_up(i,k,m,1:4)
        END DO
      END DO
    END DO

    !------------------------------------------------------------------------
    WRITE(98,'(A)') ''
    WRITE(98,'(A)') '--- OUTPUT AerZmatrices_dn (N_GEOMS x NLAYERS x 4 x 4) ---'
    DO i = 1, N_GEOMS
      DO k = 1, NLAYERS
        WRITE(98,'(A,I4,A,I4)') '  geom=', i, '  layer=', k
        DO m = 1, 4
          WRITE(98,'(A,4E16.6)') '    row: ', AerZmatrices_dn(i,k,m,1:4)
        END DO
      END DO
    END DO

    !------------------------------------------------------------------------
    WRITE(98,'(A)') ''
    WRITE(98,'(A)') '--- OUTPUT AerCoeffs (NLAYERS x 0:nmom x 6) ---'
    DO k = 1, NLAYERS
      DO m = 0, nmom
        WRITE(98,'(A,I4,A,I4,A,6E16.6)') &
          '  layer=', k, '  moment=', m, '  : ', AerCoeffs(k,m,1:6)
      END DO
    END DO

    !========================================================================
    WRITE(98,'(A)') ''
    WRITE(98,'(A)') '========================================================'
    WRITE(98,'(A)') '  END OF DUMP'
    WRITE(98,'(A)') '========================================================'

    CLOSE(98)

  end subroutine write_vfzmat_Master_debug

end module vfzmat_Master_debug_m
