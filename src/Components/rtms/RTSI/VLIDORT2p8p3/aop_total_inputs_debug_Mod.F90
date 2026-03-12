module aop_total_inputs_debug_m

  implicit none

contains

  subroutine write_aop_total_inputs_debug( &
      MAXLAYERS, MAXMOMENTS_INPUT, MAX_GEOMETRIES,  &
      NLAYERS, NSTREAMS,                             &
      N_GEOMS,                                       &
      TS_DO_SSCORR_USEFMAT,                          &
      deltau_vert_input,                             &
      omega_total_input,                             &
      greekmat_total_input,                          &
      fmatrix_up,                                    &
      fmatrix_dn )

    implicit none

    !--- Input arguments ---
    integer, intent(in) :: MAXLAYERS, MAXMOMENTS_INPUT, MAX_GEOMETRIES
    integer, intent(in) :: NLAYERS, NSTREAMS, N_GEOMS
    logical, intent(in) :: TS_DO_SSCORR_USEFMAT
    real*8,  intent(in) :: deltau_vert_input(MAXLAYERS)
    real*8,  intent(in) :: omega_total_input(MAXLAYERS)
    real*8,  intent(in) :: greekmat_total_input(0:MAXMOMENTS_INPUT, MAXLAYERS, 16)
    real*8,  intent(in) :: fmatrix_up(MAXLAYERS, MAX_GEOMETRIES, 6)
    real*8,  intent(in) :: fmatrix_dn(MAXLAYERS, MAX_GEOMETRIES, 6)

    !--- Local variables ---
    integer :: k, m, n, g

    !--- Open output file ---
    OPEN(UNIT=97, FILE='aop_total_inputs_debug.txt', STATUS='REPLACE', ACTION='WRITE')

    !========================================================================
    WRITE(97,'(A)') '========================================================'
    WRITE(97,'(A)') '   AOP POINTER ASSIGNMENT: DIAGNOSTIC DUMP'
    WRITE(97,'(A)') '========================================================'

    !------------------------------------------------------------------------
    WRITE(97,'(A)') ''
    WRITE(97,'(A)') '--- DIMENSIONS ---'
    WRITE(97,'(A,I0)') '  MAXLAYERS            = ', MAXLAYERS
    WRITE(97,'(A,I0)') '  MAXMOMENTS_INPUT     = ', MAXMOMENTS_INPUT
    WRITE(97,'(A,I0)') '  MAX_GEOMETRIES       = ', MAX_GEOMETRIES
    WRITE(97,'(A,I0)') '  NLAYERS              = ', NLAYERS
    WRITE(97,'(A,I0)') '  NSTREAMS             = ', NSTREAMS
    WRITE(97,'(A,I0)') '  N_GEOMS              = ', N_GEOMS

    !------------------------------------------------------------------------
    WRITE(97,'(A)') ''
    WRITE(97,'(A)') '--- FLAG ---'
    WRITE(97,'(A,L1)') '  TS_DO_SSCORR_USEFMAT = ', TS_DO_SSCORR_USEFMAT

    !------------------------------------------------------------------------
    WRITE(97,'(A)') ''
    WRITE(97,'(A)') '========================================================'
    WRITE(97,'(A)') '--- deltau_vert_input (MAXLAYERS) ---'
    WRITE(97,'(A)') '========================================================'
    DO k = 1, MAXLAYERS
      WRITE(97,'(A,I4,A,E20.10)') '  layer=', k, '  : ', deltau_vert_input(k)
    END DO

    !------------------------------------------------------------------------
    WRITE(97,'(A)') ''
    WRITE(97,'(A)') '========================================================'
    WRITE(97,'(A)') '--- omega_total_input (MAXLAYERS) ---'
    WRITE(97,'(A)') '========================================================'
    DO k = 1, MAXLAYERS
      WRITE(97,'(A,I4,A,E20.10)') '  layer=', k, '  : ', omega_total_input(k)
    END DO

    !------------------------------------------------------------------------
    WRITE(97,'(A)') ''
    WRITE(97,'(A)') '========================================================'
    WRITE(97,'(A)') '--- FMATRIX SECTION ---'
    WRITE(97,'(A)') '========================================================'
    WRITE(97,'(A,L1)') '  TS_DO_SSCORR_USEFMAT = ', TS_DO_SSCORR_USEFMAT

    IF ( TS_DO_SSCORR_USEFMAT ) THEN

      WRITE(97,'(A)') ''
      WRITE(97,'(A)') '  fmatrix_up (MAXLAYERS x MAX_GEOMETRIES x 6):'
      DO k = 1, MAXLAYERS
        DO g = 1, N_GEOMS
          WRITE(97,'(A,I4,A,I4,A,6E16.6)') &
            '  layer=', k, '  geom=', g, '  : ', fmatrix_up(k,g,1:6)
        END DO
      END DO

      WRITE(97,'(A)') ''
      WRITE(97,'(A)') '  fmatrix_dn (MAXLAYERS x MAX_GEOMETRIES x 6):'
      DO k = 1, MAXLAYERS
        DO g = 1, N_GEOMS
          WRITE(97,'(A,I4,A,I4,A,6E16.6)') &
            '  layer=', k, '  geom=', g, '  : ', fmatrix_dn(k,g,1:6)
        END DO
      END DO

    ELSE

      WRITE(97,'(A)') ''
      WRITE(97,'(A)') '  TS_DO_SSCORR_USEFMAT = FALSE'
      WRITE(97,'(A)') '  fmatrix_up and fmatrix_dn were NOT assigned -- skipping'

    END IF

    !------------------------------------------------------------------------
    WRITE(97,'(A)') ''
    WRITE(97,'(A)') '========================================================'
    WRITE(97,'(A)') '--- greekmat_total_input (0:MAXMOMENTS_INPUT x MAXLAYERS x 16) ---'
    WRITE(97,'(A)') '========================================================'
    DO n = 1, 16
      WRITE(97,'(A)') ''
      WRITE(97,'(A,I2,A)') '  --- Greek matrix element n=', n, ' ---'
      DO k = 1, MAXLAYERS
        WRITE(97,'(A,I4,A,I2,A)', ADVANCE='NO') '  layer=', k, '  elem=', n, '  moments: '
        DO m = 0, NSTREAMS*2+1
          WRITE(97,'(E14.6)', ADVANCE='NO') greekmat_total_input(m,k,n)
        END DO
        WRITE(97,*)
      END DO
    END DO

    !========================================================================
    WRITE(97,'(A)') ''
    WRITE(97,'(A)') '========================================================'
    WRITE(97,'(A)') '  END OF DUMP'
    WRITE(97,'(A)') '========================================================'

    CLOSE(97)

  end subroutine write_aop_total_inputs_debug

end module aop_total_inputs_debug_m
