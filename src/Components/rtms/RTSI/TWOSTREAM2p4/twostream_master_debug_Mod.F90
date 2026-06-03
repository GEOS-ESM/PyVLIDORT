module twostream_master_debug_m

  implicit none

  INTEGER, PARAMETER :: dp = KIND(1.0D0)

contains

  subroutine write_twostream_master_debug &
        ( MAXLAYERS, MAXTOTAL, MAXMESSAGES, MAXBEAMS, MAX_GEOMETRIES,     & ! Dimensions
          MAX_USER_RELAZMS, MAX_USER_STREAMS, MAX_USER_OBSGEOMS,          & ! Dimensions
          DO_UPWELLING, DO_DNWELLING, DO_PLANE_PARALLEL, DO_2S_LEVELOUT,  & ! Inputs
          DO_MVOUT_ONLY, DO_ADDITIONAL_MVOUT,                             & ! Inputs
          DO_SOLAR_SOURCES, DO_THERMAL_EMISSION, DO_SURFACE_EMISSION,     & ! Inputs
          DO_D2S_SCALING, DO_BRDF_SURFACE, DO_USER_OBSGEOMS,              & ! Inputs
          DO_SURFACE_LEAVING, DO_SL_ISOTROPIC, DO_PENTADIAG_INVERSE,      & ! Inputs
          BVPINDEX, BVPSCALEFACTOR, TAYLOR_ORDER, TAYLOR_SMALL,           & ! Inputs
          NLAYERS, NTOTAL, STREAM_VALUE, N_USER_OBSGEOMS, USER_OBSGEOMS,  & ! Inputs
          N_USER_STREAMS, USER_ANGLES, N_USER_RELAZMS, USER_RELAZMS,      & ! Inputs
          FLUX_FACTOR, NBEAMS, BEAM_SZAS, EARTH_RADIUS, HEIGHT_GRID,      & ! Inputs
          DELTAU_INPUT, OMEGA_INPUT, ASYMM_INPUT, D2S_SCALING,            & ! Inputs
          THERMAL_BB_INPUT, LAMBERTIAN_ALBEDO, BRDF_F_0, BRDF_F, UBRDF_F, & ! Inputs
          EMISSIVITY, SURFBB, SLTERM_ISOTROPIC, SLTERM_F_0 )                ! Inputs

    implicit none

    !--- 1. Dimensions ---
    integer, intent(in) :: MAXLAYERS, MAXTOTAL, MAXMESSAGES, MAXBEAMS, MAX_GEOMETRIES
    integer, intent(in) :: MAX_USER_RELAZMS, MAX_USER_STREAMS, MAX_USER_OBSGEOMS

    !--- 2. Directional & Mode Flags ---
    logical, intent(in) :: DO_UPWELLING, DO_DNWELLING, DO_PLANE_PARALLEL, DO_2S_LEVELOUT
    logical, intent(in) :: DO_MVOUT_ONLY, DO_ADDITIONAL_MVOUT
    logical, intent(in) :: DO_SOLAR_SOURCES, DO_THERMAL_EMISSION, DO_SURFACE_EMISSION
    logical, intent(in) :: DO_D2S_SCALING, DO_BRDF_SURFACE, DO_USER_OBSGEOMS
    logical, intent(in) :: DO_SURFACE_LEAVING, DO_SL_ISOTROPIC, DO_PENTADIAG_INVERSE

    !--- 3. Mathematical / Setup Scalars ---
    integer,       intent(in) :: BVPINDEX
    real(kind=dp), intent(in) :: BVPSCALEFACTOR
    integer,       intent(in) :: TAYLOR_ORDER
    real(kind=dp), intent(in) :: TAYLOR_SMALL

    !--- 4. Active Dimensions & Angles ---
    integer,       intent(in) :: NLAYERS, NTOTAL
    real(kind=dp), intent(in) :: STREAM_VALUE
    integer,       intent(in) :: N_USER_OBSGEOMS
    real(kind=dp), intent(in) :: USER_OBSGEOMS(MAX_USER_OBSGEOMS, 3)

    integer,       intent(in) :: N_USER_STREAMS
    real(kind=dp), intent(in) :: USER_ANGLES(MAX_USER_STREAMS)
    integer,       intent(in) :: N_USER_RELAZMS
    real(kind=dp), intent(in) :: USER_RELAZMS(MAX_USER_RELAZMS)

    real(kind=dp), intent(in) :: FLUX_FACTOR
    integer,       intent(in) :: NBEAMS
    real(kind=dp), intent(in) :: BEAM_SZAS(MAXBEAMS)
    
    real(kind=dp), intent(in) :: EARTH_RADIUS
    real(kind=dp), intent(in) :: HEIGHT_GRID(0:MAXLAYERS)

    !--- 5. Atmospheric Optical Properties ---
    real(kind=dp), intent(in) :: DELTAU_INPUT(MAXLAYERS)
    real(kind=dp), intent(in) :: OMEGA_INPUT(MAXLAYERS)
    real(kind=dp), intent(in) :: ASYMM_INPUT(MAXLAYERS)
    real(kind=dp), intent(in) :: D2S_SCALING(MAXLAYERS)
    real(kind=dp), intent(in) :: THERMAL_BB_INPUT(0:MAXLAYERS)

    !--- 6. Surface Properties & BRDF ---
    real(kind=dp), intent(in) :: LAMBERTIAN_ALBEDO
    real(kind=dp), intent(in) :: BRDF_F_0(0:1, MAXBEAMS)
    real(kind=dp), intent(in) :: BRDF_F(0:1)
    real(kind=dp), intent(in) :: UBRDF_F(0:1, MAX_USER_STREAMS)
    
    real(kind=dp), intent(in) :: EMISSIVITY
    real(kind=dp), intent(in) :: SURFBB
    real(kind=dp), intent(in) :: SLTERM_ISOTROPIC(MAXBEAMS)
    real(kind=dp), intent(in) :: SLTERM_F_0(0:1, MAXBEAMS)


    !--- Local variables ---
    integer :: k, b, s, u
    write(*,*) 'writing TWOSTREAM debug'
    !--- Open output file ---
    OPEN(UNIT=98, FILE='twostream_master_inputs_debug.txt', STATUS='REPLACE', ACTION='WRITE')

    !========================================================================
    WRITE(98,'(A)') '========================================================'
    WRITE(98,'(A)') '   TWOSTREAM MASTER INPUTS: DIAGNOSTIC DUMP'
    WRITE(98,'(A)') '========================================================'

    !------------------------------------------------------------------------
    WRITE(98,'(A)') ''
    WRITE(98,'(A)') '--- DIMENSIONS & ACTIVE COUNTS ---'
    WRITE(98,'(A,I0)') '  MAXLAYERS            = ', MAXLAYERS
    WRITE(98,'(A,I0)') '  MAXTOTAL             = ', MAXTOTAL
    WRITE(98,'(A,I0)') '  MAXMESSAGES          = ', MAXMESSAGES
    WRITE(98,'(A,I0)') '  MAXBEAMS             = ', MAXBEAMS
    WRITE(98,'(A,I0)') '  MAX_GEOMETRIES       = ', MAX_GEOMETRIES
    WRITE(98,'(A,I0)') '  MAX_USER_RELAZMS     = ', MAX_USER_RELAZMS
    WRITE(98,'(A,I0)') '  MAX_USER_STREAMS     = ', MAX_USER_STREAMS
    WRITE(98,'(A,I0)') '  MAX_USER_OBSGEOMS    = ', MAX_USER_OBSGEOMS
    WRITE(98,'(A)') '----------------------------------'
    WRITE(98,'(A,I0)') '  NLAYERS              = ', NLAYERS
    WRITE(98,'(A,I0)') '  NTOTAL               = ', NTOTAL
    WRITE(98,'(A,I0)') '  NBEAMS               = ', NBEAMS
    WRITE(98,'(A,I0)') '  N_USER_OBSGEOMS      = ', N_USER_OBSGEOMS
    WRITE(98,'(A,I0)') '  N_USER_STREAMS       = ', N_USER_STREAMS
    WRITE(98,'(A,I0)') '  N_USER_RELAZMS       = ', N_USER_RELAZMS

    WRITE(98,'(A)') ''
    WRITE(98,'(A)') '--- SCALARS ---'
    WRITE(98,'(A,E15.6)') '  STREAM_VALUE         = ', STREAM_VALUE
    WRITE(98,'(A,E15.6)') '  FLUX_FACTOR          = ', FLUX_FACTOR
    WRITE(98,'(A,E15.6)') '  EARTH_RADIUS         = ', EARTH_RADIUS
    WRITE(98,'(A,E15.6)') '  LAMBERTIAN_ALBEDO    = ', LAMBERTIAN_ALBEDO
    WRITE(98,'(A,E15.6)') '  SURFBB               = ', SURFBB
    WRITE(98,'(A,E15.6)') '  EMISSIVITY           = ', EMISSIVITY
    WRITE(98,'(A,I0)')    '  BVPINDEX             = ', BVPINDEX
    WRITE(98,'(A,E15.6)') '  BVPSCALEFACTOR       = ', BVPSCALEFACTOR
    WRITE(98,'(A,I0)')    '  TAYLOR_ORDER         = ', TAYLOR_ORDER
    WRITE(98,'(A,E15.6)') '  TAYLOR_SMALL         = ', TAYLOR_SMALL

    WRITE(98,'(A)') ''
    WRITE(98,'(A)') '--- LOGICAL CONTROLS ---'
    WRITE(98,'(A,L1)') '  DO_UPWELLING         = ', DO_UPWELLING
    WRITE(98,'(A,L1)') '  DO_DNWELLING         = ', DO_DNWELLING
    WRITE(98,'(A,L1)') '  DO_PLANE_PARALLEL    = ', DO_PLANE_PARALLEL
    WRITE(98,'(A,L1)') '  DO_2S_LEVELOUT       = ', DO_2S_LEVELOUT
    WRITE(98,'(A,L1)') '  DO_MVOUT_ONLY        = ', DO_MVOUT_ONLY
    WRITE(98,'(A,L1)') '  DO_ADDITIONAL_MVOUT  = ', DO_ADDITIONAL_MVOUT
    WRITE(98,'(A,L1)') '  DO_SOLAR_SOURCES     = ', DO_SOLAR_SOURCES
    WRITE(98,'(A,L1)') '  DO_THERMAL_EMISSION  = ', DO_THERMAL_EMISSION
    WRITE(98,'(A,L1)') '  DO_SURFACE_EMISSION  = ', DO_SURFACE_EMISSION
    WRITE(98,'(A,L1)') '  DO_D2S_SCALING       = ', DO_D2S_SCALING
    WRITE(98,'(A,L1)') '  DO_BRDF_SURFACE      = ', DO_BRDF_SURFACE
    WRITE(98,'(A,L1)') '  DO_USER_OBSGEOMS     = ', DO_USER_OBSGEOMS
    WRITE(98,'(A,L1)') '  DO_SURFACE_LEAVING   = ', DO_SURFACE_LEAVING
    WRITE(98,'(A,L1)') '  DO_SL_ISOTROPIC      = ', DO_SL_ISOTROPIC
    WRITE(98,'(A,L1)') '  DO_PENTADIAG_INVERSE = ', DO_PENTADIAG_INVERSE

    WRITE(98,'(A)') ''
    WRITE(98,'(A)') '--- GEOMETRY ARRAYS ---'
    if (NBEAMS > 0) then
        WRITE(98,'(A)') '  BEAM_SZAS:'
        do b = 1, NBEAMS
            WRITE(98,'(I4, E15.6)') b, BEAM_SZAS(b)
        end do
    end if

    if (N_USER_STREAMS > 0) then
        WRITE(98,'(A)') '  USER_ANGLES:'
        do s = 1, N_USER_STREAMS
            WRITE(98,'(I4, E15.6)') s, USER_ANGLES(s)
        end do
    end if

    if (N_USER_RELAZMS > 0) then
        WRITE(98,'(A)') '  USER_RELAZMS:'
        do s = 1, N_USER_RELAZMS
            WRITE(98,'(I4, E15.6)') s, USER_RELAZMS(s)
        end do
    end if

    if (N_USER_OBSGEOMS > 0) then
        WRITE(98,'(A)') '  USER_OBSGEOMS (SZA, SENSOR_ZEN, REL_AZM):'
        do u = 1, N_USER_OBSGEOMS
            WRITE(98,'(I4, 3E15.6)') u, USER_OBSGEOMS(u,1), USER_OBSGEOMS(u,2), USER_OBSGEOMS(u,3)
        end do
    end if

    WRITE(98,'(A)') ''
    WRITE(98,'(A)') '--- BRDF & SURFACE ARRAYS ---'
    WRITE(98,'(A,2E15.6)') '  BRDF_F(0:1) = ', BRDF_F(0), BRDF_F(1)
    
    if (NBEAMS > 0) then
        WRITE(98,'(A)') '  BRDF_F_0(0:1, b):'
        do b = 1, NBEAMS
            WRITE(98,'(I4, 2E15.6)') b, BRDF_F_0(0,b), BRDF_F_0(1,b)
        end do
        WRITE(98,'(A)') '  SLTERM_ISOTROPIC(b):'
        do b = 1, NBEAMS
            WRITE(98,'(I4, E15.6)') b, SLTERM_ISOTROPIC(b)
        end do
        WRITE(98,'(A)') '  SLTERM_F_0(0:1, b):'
        do b = 1, NBEAMS
            WRITE(98,'(I4, 2E15.6)') b, SLTERM_F_0(0,b), SLTERM_F_0(1,b)
        end do
    end if

    if (N_USER_STREAMS > 0) then
        WRITE(98,'(A)') '  UBRDF_F(0:1, s):'
        do s = 1, N_USER_STREAMS
            WRITE(98,'(I4, 2E15.6)') s, UBRDF_F(0,s), UBRDF_F(1,s)
        end do
    end if

    WRITE(98,'(A)') ''
    WRITE(98,'(A)') '--- LAYER INPUT ARRAYS (1 to NLAYERS) ---'
    WRITE(98,'(A)') '  k    HEIGHT_GRID   DELTAU_INPUT  OMEGA_INPUT   ASYMM_INPUT   D2S_SCALING   THERMAL_BB'
    
    do k = 1, NLAYERS
      WRITE(98,'(I4, 6E14.5)') k, HEIGHT_GRID(k), DELTAU_INPUT(k), &
                               OMEGA_INPUT(k), ASYMM_INPUT(k), &
                               D2S_SCALING(k), THERMAL_BB_INPUT(k)
    end do
    
    ! Write out level 0 (Surface) for 0-indexed arrays
    WRITE(98,'(A,2E14.5)') ' Sfc 0 HEIGHT_GRID & THERMAL_BB = ', HEIGHT_GRID(0), THERMAL_BB_INPUT(0)

    WRITE(98,'(A)') '========================================================'
    CLOSE(98)

  end subroutine write_twostream_master_debug

end module twostream_master_debug_m
