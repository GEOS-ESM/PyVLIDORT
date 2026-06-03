module twostream_master_outputs_debug_m

  implicit none

  INTEGER, PARAMETER :: dp = KIND(1.0D0)

contains

  subroutine write_twostream_master_outputs_debug &
        ( MAXLAYERS, MAXMESSAGES, MAXBEAMS, MAX_GEOMETRIES,       & ! Dimensions needed for array bounds
          NLAYERS, NBEAMS, N_GEOMETRIES,                          & ! Active dimensions for looping
          INTENSITY_TOA, INTENSITY_BOA, FLUXES_TOA, FLUXES_BOA,   & ! Radiance & Flux Outputs
          RADLEVEL_UP, RADLEVEL_DN,                               & ! Level Output
          STATUS_INPUTCHECK, C_NMESSAGES, C_MESSAGES, C_ACTIONS,  & ! Exception handling
          STATUS_EXECUTION, E_MESSAGE, E_TRACE_1, E_TRACE_2 )       ! Exception handling

    implicit none

    !--- Dimensions ---
    integer, intent(in) :: MAXLAYERS, MAXMESSAGES, MAXBEAMS, MAX_GEOMETRIES
    integer, intent(in) :: NLAYERS, NBEAMS, N_GEOMETRIES

    !--- Radiance Results ---
    real(kind=dp), intent(in) :: INTENSITY_TOA(MAX_GEOMETRIES)
    real(kind=dp), intent(in) :: INTENSITY_BOA(MAX_GEOMETRIES)

    !--- Flux Output (MAXBEAMS, 2: usually upwelling and downwelling/direct) ---
    real(kind=dp), intent(in) :: FLUXES_TOA(MAXBEAMS, 2)
    real(kind=dp), intent(in) :: FLUXES_BOA(MAXBEAMS, 2)

    !--- Optional Output at ALL LEVELS ---
    real(kind=dp), intent(in) :: RADLEVEL_UP(MAX_GEOMETRIES, 0:MAXLAYERS)
    real(kind=dp), intent(in) :: RADLEVEL_DN(MAX_GEOMETRIES, 0:MAXLAYERS)

    !--- Exception handling ---
    integer,       intent(in) :: STATUS_INPUTCHECK, C_NMESSAGES, STATUS_EXECUTION
    character*100, intent(in) :: C_MESSAGES(0:MAXMESSAGES)
    character*100, intent(in) :: C_ACTIONS(0:MAXMESSAGES)
    character*100, intent(in) :: E_MESSAGE, E_TRACE_1, E_TRACE_2

    !--- Local variables ---
    integer :: g, b, k, msg

    !--- Open output file ---
    OPEN(UNIT=99, FILE='twostream_master_outputs_debug.txt', STATUS='REPLACE', ACTION='WRITE')

    !========================================================================
    WRITE(99,'(A)') '========================================================'
    WRITE(99,'(A)') '   TWOSTREAM MASTER OUTPUTS: DIAGNOSTIC DUMP'
    WRITE(99,'(A)') '========================================================'

    !------------------------------------------------------------------------
    WRITE(99,'(A)') ''
    WRITE(99,'(A)') '--- ACTIVE DIMENSIONS ---'
    WRITE(99,'(A,I0)') '  N_GEOMETRIES (Returned) = ', N_GEOMETRIES
    WRITE(99,'(A,I0)') '  NBEAMS                  = ', NBEAMS
    WRITE(99,'(A,I0)') '  NLAYERS                 = ', NLAYERS

    WRITE(99,'(A)') ''
    WRITE(99,'(A)') '--- EXCEPTION HANDLING ---'
    WRITE(99,'(A,I0)') '  STATUS_INPUTCHECK       = ', STATUS_INPUTCHECK
    WRITE(99,'(A,I0)') '  C_NMESSAGES             = ', C_NMESSAGES
    
    if (C_NMESSAGES > 0) then
      WRITE(99,'(A)') '  Input Check Messages:'
      do msg = 1, min(C_NMESSAGES, MAXMESSAGES)
        WRITE(99,'(A,I0,A,A)') '    Msg ', msg, ': ', TRIM(C_MESSAGES(msg))
        WRITE(99,'(A,I0,A,A)') '    Act ', msg, ': ', TRIM(C_ACTIONS(msg))
      end do
    end if

    WRITE(99,'(A,I0)') '  STATUS_EXECUTION        = ', STATUS_EXECUTION
    if (STATUS_EXECUTION /= 0) then
      WRITE(99,'(A,A)') '  E_MESSAGE               = ', TRIM(E_MESSAGE)
      WRITE(99,'(A,A)') '  E_TRACE_1               = ', TRIM(E_TRACE_1)
      WRITE(99,'(A,A)') '  E_TRACE_2               = ', TRIM(E_TRACE_2)
    end if

    WRITE(99,'(A)') ''
    WRITE(99,'(A)') '--- INTENSITIES (TOA & BOA) ---'
    if (N_GEOMETRIES > 0) then
      WRITE(99,'(A)') '  Geom    INTENSITY_TOA      INTENSITY_BOA'
      do g = 1, N_GEOMETRIES
        WRITE(99,'(I6, 2E19.8)') g, INTENSITY_TOA(g), INTENSITY_BOA(g)
      end do
    else
      WRITE(99,'(A)') '  (N_GEOMETRIES is 0 - no intensities to display)'
    end if

    WRITE(99,'(A)') ''
    WRITE(99,'(A)') '--- FLUXES (TOA & BOA) ---'
    if (NBEAMS > 0) then
      WRITE(99,'(A)') '  Beam    FLUX_TOA(1)        FLUX_TOA(2)        FLUX_BOA(1)        FLUX_BOA(2)'
      do b = 1, NBEAMS
        WRITE(99,'(I6, 4E19.8)') b, FLUXES_TOA(b,1), FLUXES_TOA(b,2), FLUXES_BOA(b,1), FLUXES_BOA(b,2)
      end do
    else
      WRITE(99,'(A)') '  (NBEAMS is 0 - no fluxes to display)'
    end if

    WRITE(99,'(A)') ''
    WRITE(99,'(A)') '--- RADLEVELS (ALL LEVELS) ---'
    if (N_GEOMETRIES > 0 .and. NLAYERS > 0) then
      ! Loop over geometries first, then levels, to make it readable
      do g = 1, N_GEOMETRIES
        WRITE(99,'(A,I0,A)') '  Geometry ', g, ':'
        WRITE(99,'(A)')      '    Level    RADLEVEL_UP        RADLEVEL_DN'
        ! Level 0 (Surface / TOA depending on convention)
        WRITE(99,'(4X, I5, 2E19.8)') 0, RADLEVEL_UP(g, 0), RADLEVEL_DN(g, 0)
        do k = 1, NLAYERS
          WRITE(99,'(4X, I5, 2E19.8)') k, RADLEVEL_UP(g, k), RADLEVEL_DN(g, k)
        end do
        WRITE(99,'(A)') ''
      end do
    else
      WRITE(99,'(A)') '  (N_GEOMETRIES or NLAYERS is 0 - no radlevels to display)'
    end if

    WRITE(99,'(A)') '========================================================'
    CLOSE(99)

  end subroutine write_twostream_master_outputs_debug

end module twostream_master_outputs_debug_m
