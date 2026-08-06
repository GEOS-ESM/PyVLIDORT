PROGRAM test_2s_ONLY_2p4

  USE twostream_brdf_supplement_m
  USE twostream_master_m
  USE TWOSTREAM_getPlanck

      implicit none

!  Notes 21 december 2012. Observational Geometry Inputs. Marked with !@@

!     Observation-Geometry New dimensioning.    MAX_USER_OBSGEOMS
!     Observation-Geometry input control.       DO_USER_OBSGEOMS
!     Observation-Geometry input control.       N_USER_OBSGEOMS
!     User-defined Observation Geometry angles. USER_OBSGEOMS

!  17 july    2013, ALL LEVEL output,      Version 2.2
!  December   2013, Flux outputs + control Version 2.3
!  23 january 2014, Surface-leaving        Version 2.3
!  25 June    2014, BVProblem control      Version 2.3
!  15 August  2014, Beam/Thermal Greens    Version 2.4
!  08 January 2015, Greens, Linearized     Version 2.4

!  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!          2STREAM ARGUMENTS
!  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  precision

      INTEGER, PARAMETER :: dp = KIND( 1.0D0 )

!  Dimensioning integer parameters

      INTEGER, PARAMETER :: max_user_obsgeoms = 2     !@@
      INTEGER, PARAMETER :: maxbeams         = max_user_obsgeoms
      INTEGER, PARAMETER :: max_user_angles  = max_user_obsgeoms
      INTEGER, PARAMETER :: max_user_relazms = max_user_obsgeoms

      INTEGER, PARAMETER :: maxmessages      = 25
      INTEGER, PARAMETER :: maxlayers        = 23

      INTEGER, PARAMETER :: maxstreams_brdf     = 50
      INTEGER, PARAMETER :: max_brdf_kernels    = 3
      INTEGER, PARAMETER :: max_brdf_parameters = 3

      INTEGER, PARAMETER :: max_geometries = maxbeams * max_user_angles * max_user_relazms
      INTEGER, PARAMETER :: maxtotal       = 2*maxlayers

!  Directional Flags

      LOGICAL       :: DO_UPWELLING, DO_DNWELLING

!  Plane parallel and deltam-2stream scaling flags

      LOGICAL       :: DO_PLANE_PARALLEL, DO_D2S_SCALING

!  BRDF surface flag

      LOGICAL       :: DO_BRDF_SURFACE

!   !@@ Observational Geometry flag !@@ 2p1

      LOGICAL       :: DO_USER_OBSGEOMS !@@ 2p1

!  @@ Rob Spurr, 17 July 2013, Version 2.2, Levelout flag

      LOGICAL       :: DO_2S_LEVELOUT     ! @@ 2p2

!  @@ Rob Spurr, 05 November 2013, Version 2.3, Flux option flags

      LOGICAL       :: DO_MVOUT_ONLY       ! @@ 2p3
      LOGICAL       :: DO_ADDITIONAL_MVOUT ! @@ 2p3

!  @@ Rob Spurr, 23 January 2014, Version 2.3, SLEAVE option flags

      LOGICAL       :: DO_SURFACE_LEAVING
      LOGICAL       :: DO_SL_ISOTROPIC

!  ** New **. October 2011, Sources control, including thermal

      LOGICAL       :: DO_THERMAL_EMISSION
      LOGICAL       :: DO_SURFACE_EMISSION
      LOGICAL       :: DO_SOLAR_SOURCES

!  Order of Taylor series (including terms up to EPS^n). @@ 2p4
      
      INTEGER       :: TAYLOR_ORDER
      REAL(kind=dp) :: TAYLOR_SMALL

!  BVP control --- New 6/25/14, Version 2.3 and higher
!   * PentaDiagonal Inverse flag (BVP solved from bottom to top). Only for BVPIndex = 1
!   * BVP Index : 0 = LAPACK, 1 = Penta # 1 (original), 2 = Penta # 2 (new, 2012 Kanal paper)
!   * BVP Scale Factor. Debug only. Set this to 1.0 on input

      LOGICAL        :: DO_PENTADIAG_INVERSE
      INTEGER        :: BVPINDEX
      REAL(kind=dp)  :: BVPSCALEFACTOR

!  Numbers

      INTEGER       :: NLAYERS, NTOTAL
      INTEGER       :: N_USER_ANGLES, N_USER_RELAZMS, NBEAMS

!  Geometry

      REAL(KIND=DP) :: BEAM_SZAS    ( MAXBEAMS )
      REAL(KIND=DP) :: USER_ANGLES  ( MAX_USER_ANGLES )
      REAL(KIND=DP) :: USER_RELAZMS ( MAX_USER_RELAZMS )

!  Stream value

      REAL(KIND=DP) :: STREAM_VALUE

!  Observational geometry input. [Same as LIDORT]. New 12/21/12 !@@

      INTEGER       :: N_USER_OBSGEOMS                    !@@
      REAL(kind=dp) :: USER_OBSGEOMS(MAX_USER_OBSGEOMS,3) !@@

!  Thermal inputs

      REAL(KIND=DP) :: SURFBB
      REAL(KIND=DP) :: THERMAL_BB_INPUT ( 0:MAXLAYERS )

!  Lambertian Surface control (threaded)

      REAL(KIND=DP) :: LAMBERTIAN_ALBEDO

!  BRDF Fourier components (NOT threaded)
!  0 and 1 Fourier components of BRDF, following order (same all threads)
!    incident solar directions,  reflected quadrature stream
!    incident quadrature stream, reflected quadrature stream
!    incident solar directions,  reflected user streams
!    incident quadrature stream, reflected user streams

      REAL(KIND=DP) :: BRDF_F_0 ( 0:1, MAXBEAMS )
      REAL(KIND=DP) :: BRDF_F   ( 0:1 )
!     REAL(KIND=DP) :: UBRDF_F_0 ( 0:1, MAX_USER_ANGLES, MAXBEAMS )
      REAL(KIND=DP) :: UBRDF_F  ( 0:1, MAX_USER_ANGLES )

!  Emissivity

      REAL(KIND=DP) :: EMISSIVITY

!  Version 2p3. 1/23/14. Introduce SLEAVE stuff
!    Do not require any first-order inputs (exact or Fourier)
!    Isotropic Surface leaving term (if flag set)
!  Fourier components of Surface-leaving terms:

      REAL(kind=dp) ::  SLTERM_ISOTROPIC ( MAXBEAMS )
      REAL(kind=dp) ::  SLTERM_F_0 ( 0:1, MAXBEAMS )

!  Flux factor

      REAL(KIND=DP) :: FLUX_FACTOR

!  height and earth radius

      REAL(KIND=DP) :: EARTH_RADIUS
      REAL(KIND=DP) :: HEIGHT_GRID ( 0:MAXLAYERS )

!  Atmospheric Optical properties

      REAL(KIND=DP) :: DELTAU_INPUT( MAXLAYERS )
      REAL(KIND=DP) :: OMEGA_INPUT ( MAXLAYERS )
      REAL(KIND=DP) :: ASYMM_INPUT ( MAXLAYERS )
      REAL(KIND=DP) :: D2S_SCALING ( MAXLAYERS )

!  Results

      REAL(KIND=DP) :: INTENSITY_TOA ( MAX_GEOMETRIES )
      REAL(KIND=DP) :: INTENSITY_BOA ( MAX_GEOMETRIES )

!  Flux output
!     ! @@ Rob Spurr, 05 November 2013, Version 2.3 --> Flux Output

     REAL(kind=dp) :: FLUXES_TOA(MAXBEAMS,2)
     REAL(kind=dp) :: FLUXES_BOA(MAXBEAMS,2)

!  output solutions at ALL levels
!     ! @@ Rob Spurr, 17 July 2013, Version 2.2 --> Optional Output at ALL LEVELS

      REAL(kind=dp) :: RADLEVEL_UP (MAX_GEOMETRIES,0:MAXLAYERS)
      REAL(kind=dp) :: RADLEVEL_DN (MAX_GEOMETRIES,0:MAXLAYERS)

!  Numbers (geometry)
!   N_GEOMETRIES = NBEAMS * N_USER_STREAMS * N_USER_RELAZMS (Lattice value)
!   N_GEOMETRIES = N_USER_OBSGEOMS                          (OBsGeom value)

      INTEGER       :: n_geometries

!  Exception handling

!    1. Up to 100 Check Messages and actions

      INTEGER       :: STATUS_INPUTCHECK
      INTEGER       :: C_NMESSAGES
      CHARACTER (LEN=100) :: C_MESSAGES ( 0:MAXMESSAGES )
      CHARACTER (LEN=100) :: C_ACTIONS  ( 0:MAXMESSAGES )

!    2. Execution message and 2 Traces

      INTEGER       :: STATUS_EXECUTION
      CHARACTER (LEN=100) :: E_MESSAGE, E_TRACE_1, E_TRACE_2

!  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!          OTHER ARGUMENTS
!  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  External thread parameter

      INTEGER, PARAMETER :: maxthreads = 6

!  Saved values

      REAL(KIND=dp) :: INTENSITY_TOA_Save (MAX_GEOMETRIES,MAXTHREADS)
      REAL(KIND=dp) :: INTENSITY_BOA_Save (MAX_GEOMETRIES,MAXTHREADS)
      REAL(kind=dp) :: FLUXES_TOA_Save    (MAXBEAMS,2,MAXTHREADS)
      REAL(kind=dp) :: FLUXES_BOA_Save    (MAXBEAMS,2,MAXTHREADS)
      REAL(kind=dp) :: RADLEVEL_UP_Save   (MAX_GEOMETRIES,0:MAXLAYERS,MAXTHREADS)
      REAL(kind=dp) :: RADLEVEL_DN_Save   (MAX_GEOMETRIES,0:MAXLAYERS,MAXTHREADS)

!  Other flags

      LOGICAL       :: DO_FULLQUADRATURE

!  Output strings

      CHARACTER (LEN=10) :: c7
      CHARACTER (LEN=60) :: cstream

!  Help variables

      INTEGER       :: n,n6,ndum,ldum,t,k,g,v,nthreads,thread
      REAL(KIND=DP) :: kd, gaer, waer, taer, parcel, raywt, aerwt
      REAL(KIND=DP) :: aersca, aerext, molsca, totsca, totext
      REAL(KIND=DP) :: molomg ( maxlayers ), molext ( maxlayers ), m1, m2
      REAL(KIND=DP) :: raymoms ( 0:2,maxlayers )
      REAL(KIND=DP) :: albedo_save
      data albedo_save / 0.2d0 /
      REAL(KIND=DP) :: windspeed_save
      data windspeed_save / 10.d0 /
      LOGICAL       :: do_old_outformat

!  Timing variables

      logical, parameter :: do_timing = .false.
      real    :: e1, e2, exectime(maxthreads),execall(20,6,maxthreads)
      integer :: k0, m, L, nruns(6),irun,ND,LD,trun

!  Thermal emission setups

      LOGICAL       :: THERMFAIL
      CHARACTER (LEN=70) :: THERMALMESSAGE
      REAL(KIND=DP) :: SURFTEMP, AIRTEMPS ( 0:MAXLAYERS )
      REAL(KIND=DP) :: WNUMLO, WNUMHI
      INTEGER       :: NTEMPS, SMALLV

!  Cox-Munk Surface control (threaded)

      REAL(KIND=DP) :: WINDSPEED ( MAXTHREADS )

!  BRDF inputs

      INTEGER       :: NSTREAMS_BRDF
      LOGICAL       :: LAMBERTIAN_KERNEL_FLAG ( MAX_BRDF_KERNELS )
      LOGICAL       :: DO_SHADOW_EFFECT
      INTEGER       :: N_BRDF_KERNELS
      INTEGER       :: WHICH_BRDF   ( MAX_BRDF_KERNELS )
      REAL(KIND=DP) :: BRDF_FACTORS ( MAX_BRDF_KERNELS )
      INTEGER       :: N_BRDF_PARAMETERS ( MAX_BRDF_KERNELS )
      REAL(KIND=DP) :: BRDF_PARAMETERS   ( MAX_BRDF_KERNELS, MAX_BRDF_PARAMETERS )

!  BRDF Exception handling. !@@ Added, 12/31/12

      INTEGER       :: status_brdfsup
      CHARACTER*100 :: message_brdf, action_brdf

!  Start program

!  Set up inputs

      DO_UPWELLING      = .true.
      DO_DNWELLING      = .true.

!  Control read

      open(44,file = '2s_2p4_driver.inp',status = 'old')
      read(44,*)DO_SOLAR_SOURCES
      read(44,*)DO_THERMAL_EMISSION
      read(44,*)DO_SURFACE_EMISSION
      read(44,*)DO_BRDF_SURFACE
      read(44,*)DO_FULLQUADRATURE
      read(44,*)DO_D2S_SCALING
      read(44,*)DO_PLANE_PARALLEL
      read(44,*)DO_USER_OBSGEOMS ! @@ 2p1
      read(44,*)
      read(44,*)
      read(44,*)DO_2S_LEVELOUT       ! @@ 2p2
      read(44,*)DO_MVOUT_ONLY        ! @@ 2p3
      read(44,*)DO_ADDITIONAL_MVOUT  ! @@ 2p3
      read(44,*)DO_PENTADIAG_INVERSE ! @@ 2p3
      close(44)

!  Modify test file input

!      DO_PLANE_PARALLEL = .false.
      DO_BRDF_SURFACE   = .true.

!  Define other control inputs not in the input file
!  (Note: surface-leaving not active yet)

      DO_SURFACE_LEAVING = .false.
      DO_SL_ISOTROPIC    = .false.
      SLTERM_ISOTROPIC   = 0.0D0
      SLTERM_F_0         = 0.0D0

      !For SLEAVE input quick check only
      !DO_SURFACE_LEAVING = .true.
      !DO_SL_ISOTROPIC    = .true.
      !SLTERM_ISOTROPIC   = 1.0D-1
      !SLTERM_F_0         = 1.0D-1

!  Taylor control

      TAYLOR_SMALL = 1.0d-03
      TAYLOR_ORDER = 3

!  Scale factor and BVP Index are pre-set

      BVPINDEX       = 1
      BVPSCALEFACTOR = 1.0d0

!  Control integers

      nthreads = 6

      n_user_obsgeoms = max_user_obsgeoms
      nbeams          = n_user_obsgeoms
      n_user_angles   = n_user_obsgeoms
      n_user_relazms  = n_user_obsgeoms

      nlayers  = 23
      ntotal   = 2 * nlayers

      nstreams_brdf  = 50

!  Flux factor

      flux_factor = 1.0d0

!  Angles

      beam_szas(1)   = 35.0d0
      user_angles(1) = 30.0d0
      user_relazms(1)= 10.0d0
      beam_szas(2)   = 45.0d0
      user_angles(2) = 40.0d0
      user_relazms(2)= 70.0d0
!      beam_szas(3)   = 65.0d0
!      user_angles(3) = 50.0d0
!      user_relazms(3)= 170.0d0
!      beam_szas(4)   = 66.0d0
!      user_angles(4) = 51.0d0
!      user_relazms(4)= 171.0d0
!      beam_szas(5)   = 67.0d0
!      user_angles(5) = 52.0d0
!      user_relazms(5)= 172.0d0
!      beam_szas(6)   = 67.1d0
!      user_angles(6) = 52.1d0
!      user_relazms(6)= 172.1d0

      do n = 1, n_user_obsgeoms
         user_obsgeoms(n,1) = beam_szas(n)
         user_obsgeoms(n,2) = user_angles(n)
         user_obsgeoms(n,3) = user_relazms(n)
      enddo

!  Atmosphere

      earth_radius = 6371.0d0
      height_grid = 0.0d0
      open(45,file='data/input_atmos.dat',status='old' )
      read(45,'(i5,1p25e18.9)')ldum, (raymoms(0,n),n=1,nlayers)
      read(45,'(i5,1p25e18.9)')ldum, (raymoms(1,n),n=1,nlayers)
      read(45,'(i5,1p25e18.9)')ldum, (raymoms(2,n),n=1,nlayers)
      height_grid(0) = 60.0d0
      do n = 1, nlayers
         read(45,'(i4,f12.5,1p6e16.7)')ndum,height_grid(n), &
            molext(n),molomg(n),kd,kd,kd,kd
      enddo
      close(45)

!  Add Aerosols bottom 6 layers, spread evenly

      n6 = nlayers - 6; gaer = 0.8d0 ; waer = 0.95d0 ; taer = 0.5d0

!  Thermal stuff
!  -------------

!   Array of temperatures should be 24 levels

      open(1,file= 'data/input_temp23.dat', status='old')
      read(1,*)SURFTEMP
      read(1,*)NTEMPS
      do n = 0, NTEMPS
         read(1,*)ndum, AIRTEMPS(NTEMPS-N)
      enddo
      close(1)

!  Zeroing

      lambertian_albedo = 0.0d0

      deltau_input = 0.0d0
      omega_input = 0.0d0
      asymm_input = 0.0d0
      BRDF_F_0  = 0.0d0
      BRDF_F  = 0.0d0
      UBRDF_F = 0.0d0

!  Two choices of stream value................ CHOOSE One !!!!!

      if ( do_fullquadrature ) then
        stream_value = dsqrt(1.0d0 / 3.0d0)
      else
        stream_value = 0.5d0
      endif

!  Headers

      if ( stream_value .eq. 0.5d0 ) then
        c7 = 'SV1_QS_XXX'
        if (do_d2s_scaling) c7 = 'SV1_QS_D2S'
        if (do_d2s_scaling) cstream =  &
            '(Stream value = 0.50000), with delta-2 stream scaling'
        if (.not.do_d2s_scaling) cstream =  &
            '(Stream value = 0.50000), with NO delta-2 stream scaling'
      else
        c7 = 'SV2_QS_XXX'
        if (do_d2s_scaling) c7 = 'SV2_QS_D2S'
        cstream = '(Stream value = 0.57735)'
        if (do_d2s_scaling) cstream = &
            '(Stream value = 0.57735), with delta-2 stream scaling'
        if (.not.do_d2s_scaling) cstream = &
            '(Stream value = 0.57735), with NO delta-2 stream scaling'
      endif
      if ( do_plane_parallel)c7(5:6) = 'PP'

!     Initialize output

      INTENSITY_TOA_SAVE = 0.0d0
      INTENSITY_BOA_SAVE = 0.0d0
      FLUXES_TOA_SAVE    = 0.0d0
      FLUXES_BOA_SAVE    = 0.0d0
      RADLEVEL_UP_SAVE   = 0.0d0
      RADLEVEL_DN_SAVE   = 0.0d0

!  Baseline calculation 1 : RADIANCE ONLY with Basic Master
!  ========================================================

!  Start thread loop

!      do thread = 1, 1
      do thread = 1, nthreads
        t = thread

!  Thread 1: 2S, Solar only   , no   delta-M scaling
!  Thread 2: 2S, Solar only   , with delta-M scaling
!  Thread 3: 2S, Thermal only , no   delta-M scaling
!  Thread 4: 2S, Thermal only , with delta-M scaling
!  Thread 5: 2S, Solar+Thermal, no   delta-M scaling
!  Thread 6: 2S, Solar+Thermal, with delta-M scaling

        if ( thread .eq. 1 ) then
          DO_SOLAR_SOURCES    = .true.
          DO_D2S_SCALING      = .false.
          DO_THERMAL_EMISSION = .false.
          DO_SURFACE_EMISSION = .false.
        else if ( thread .eq. 2 ) then
          DO_SOLAR_SOURCES    = .true.
          DO_D2S_SCALING      = .true.
          DO_THERMAL_EMISSION = .false.
          DO_SURFACE_EMISSION = .false.
        else if ( thread .eq. 3 ) then
          DO_SOLAR_SOURCES    = .false.
          DO_D2S_SCALING      = .false.
          DO_THERMAL_EMISSION = .true.
          DO_SURFACE_EMISSION = .true.
          nbeams = 1;  n_user_relazms = 1  ! Must set these 2, otherwise doesn't work
        else if ( thread .eq. 4 ) then
          DO_SOLAR_SOURCES    = .false.
          DO_D2S_SCALING      = .true.
          DO_THERMAL_EMISSION = .true.
          DO_SURFACE_EMISSION = .true.
          nbeams = 1;  n_user_relazms = 1  ! Must set these 2, otherwise doesn't work
        else if ( thread .eq. 5 ) then
          DO_SOLAR_SOURCES    = .true.
          DO_D2S_SCALING      = .false.
          DO_THERMAL_EMISSION = .true.
          DO_SURFACE_EMISSION = .true.
          nbeams = n_user_obsgeoms;  n_user_relazms = n_user_obsgeoms  ! Must set these 2, otherwise doesn't work
        else if ( thread .eq. 6 ) then
          DO_SOLAR_SOURCES    = .true.
          DO_THERMAL_EMISSION = .true.
          DO_SURFACE_EMISSION = .true.
          DO_D2S_SCALING      = .true.
          nbeams = n_user_obsgeoms;  n_user_relazms = n_user_obsgeoms  ! Must set these 2, otherwise doesn't work
        endif
        write(*,*)'Doing thread # ', thread

!  Thermal input

        wnumhi = 2500.0d0 + 1.0d0
        wnumlo = 2500.0d0 - 1.0d0

        if ( do_thermal_emission ) then
          do n = 0, nlayers
            call get_planckfunction &
             ( wnumlo, wnumhi, AIRTEMPS(N), &
               THERMAL_BB_INPUT(N), SMALLV, THERMFAIL, THERMALMESSAGE )
          enddo
          if ( THERMFAIL ) THEN
            write(*,*)THERMALMESSAGE ; stop
          endif
        else
          thermal_bb_input = 0.0d0
        endif

        if ( do_surface_emission ) then
          call get_planckfunction &
           ( wnumlo, wnumhi, SURFTEMP, &
             SURFBB, SMALLV, THERMFAIL, THERMALMESSAGE )
          if ( THERMFAIL ) THEN
            write(*,*)THERMALMESSAGE ; stop
          endif
          emissivity = 1.0d0 - albedo_save
        else
          surfbb = 0.0d0
          emissivity = 0.0d0
        end if

!  Create optical properties, every thread

        do n = 1, n6
          deltau_input(n) = molext(n)
          omega_input(n)  = molomg(n)
          asymm_input(n)  = 0.0d0
          molsca = molomg(n) * molext(n)
          m1     = raymoms(2,n) * molsca ; m2 = 0.0d0
          d2s_scaling(n)  = ( m1 + m2 ) / molsca / 5.0d0
        enddo
        parcel = taer / ( height_grid(n6) - height_grid(nlayers) )
        do n = n6 + 1, nlayers
          aerext = Parcel * ( height_grid(n-1) - height_grid(n) )
          aersca = aerext * waer
          molsca = molomg(n) * molext(n)
          totext = molext(n) + aerext
          totsca = molsca    + aersca
          raywt  = molsca / totsca
          aerwt  = aersca / totsca
          deltau_input(n) = totext
          omega_input(n)  = totsca / totext
          asymm_input(n)  = gaer * aersca / totsca 
          m1 = raymoms(2,n) * molsca
          m2 = 5.0d0 * gaer * gaer * aersca
          d2s_scaling(n)  = ( m1 + m2 ) / totsca / 5.0d0
        enddo

!  Set up surface

        if ( DO_BRDF_SURFACE ) then
          !BRDF surface

          !(1) Initialize:

          !Basic BRDF supplement inputs
          LAMBERTIAN_KERNEL_FLAG    = .FALSE.
          DO_SHADOW_EFFECT          = .FALSE.

          N_BRDF_KERNELS            = 0

          WHICH_BRDF                = 0
          BRDF_FACTORS              = 0.0D0
          N_BRDF_PARAMETERS         = 0
          BRDF_PARAMETERS           = 0.0D0

          !(2) Now, for the Cox-Munk BRDF, define:

          !Basic BRDF supplement inputs
          N_BRDF_KERNELS            = 1

          WHICH_BRDF(1)             = 9 !Cox-Munk
          BRDF_FACTORS(1)           = 1.0D0
          N_BRDF_PARAMETERS(1)      = 3
          WINDSPEED(T)              = windspeed_save
          BRDF_PARAMETERS(1,1)      = 0.003d0+0.00512d0*WINDSPEED(T) !Slope squared
          BRDF_PARAMETERS(1,2)      = 1.33D0*1.33D0 !Square of refractive index
          !Note: The Cox-Munk shadow effect parameter (the 3rd
          !      parameter BRDF_PARAMETERS(1,3)) is set internally by
          !      TWOSTREAM_BRDFMASTER using the DO_SHADOW_EFFECT flag

          !Call BRDF supplement
          CALL TWOSTREAM_BRDFMASTER &
            ( MAXBEAMS, MAX_USER_ANGLES, MAX_USER_OBSGEOMS, & ! Dimensions !@@
              MAXSTREAMS_BRDF, MAX_BRDF_KERNELS,            & ! Dimensions
              MAX_BRDF_PARAMETERS,                          & ! Dimensions
              DO_SOLAR_SOURCES, DO_USER_OBSGEOMS,           & ! Inputs !@@
              LAMBERTIAN_KERNEL_FLAG,                       & ! Inputs
              DO_SHADOW_EFFECT, DO_SURFACE_EMISSION,        & ! Inputs
              NBEAMS, N_USER_ANGLES, N_USER_OBSGEOMS,       & ! Inputs !@@
              BEAM_SZAS, USER_ANGLES, USER_OBSGEOMS,        & ! Inputs !@@
              STREAM_VALUE, NSTREAMS_BRDF,                  & ! Inputs
              N_BRDF_KERNELS, WHICH_BRDF, BRDF_FACTORS,     & ! Inputs
              N_BRDF_PARAMETERS, BRDF_PARAMETERS,           & ! Inputs
              BRDF_F_0, BRDF_F, UBRDF_F, EMISSIVITY,        & ! Outputs
              STATUS_BRDFSUP, MESSAGE_BRDF, ACTION_BRDF )     ! Outputs

        else
          !Lambertian surface

          LAMBERTIAN_ALBEDO = albedo_save
        end if

!  Exception handling on BRDF

        if ( DO_BRDF_SURFACE ) then
           if ( STATUS_BRDFSUP .ne. 0 ) then
              write(*,'(a,i4)')'BRDF supplement Check failed from Baseline Run # ',T
              write(*,*)' - Print 1 Message and 1 Action'
              write(*,'(A)') TRIM(MESSAGE_BRDF)
              write(*,'(A)') TRIM(ACTION_BRDF)
              stop'Test_2S_only program aborted'
           endif
        endif

!  initialize timing

        trun = 1 ; if ( do_timing) trun = 1999
        exectime(thread) = 0.0
        do irun = 1, trun
        if (mod(irun,2000).eq.0)write(*,*)irun
        call cpu_time(e1)

!  Call to Regular model

        CALL TWOSTREAM_MASTER &
        ( MAXLAYERS, MAXTOTAL, MAXMESSAGES, MAXBEAMS, MAX_GEOMETRIES,     & ! Dimensions
          MAX_USER_RELAZMS, MAX_USER_ANGLES, MAX_USER_OBSGEOMS,           & ! Dimensions !@@ 2p1
          DO_UPWELLING, DO_DNWELLING, DO_PLANE_PARALLEL, DO_2S_LEVELOUT,  & ! Inputs     !@@ 2p2
          DO_MVOUT_ONLY, DO_ADDITIONAL_MVOUT,                             & ! Inputs     !@@ 2p3
          DO_SOLAR_SOURCES, DO_THERMAL_EMISSION, DO_SURFACE_EMISSION,     & ! Inputs
          DO_D2S_SCALING, DO_BRDF_SURFACE, DO_USER_OBSGEOMS,              & ! Inputs     !@@ 2p1
          DO_SURFACE_LEAVING, DO_SL_ISOTROPIC, DO_PENTADIAG_INVERSE,      & ! Input !@@ 2p3 6/25/14
          BVPINDEX, BVPSCALEFACTOR, TAYLOR_ORDER, TAYLOR_SMALL,           & ! Input !@@ 2p3 6/25/14, 8/15/14
          NLAYERS, NTOTAL, STREAM_VALUE, N_USER_OBSGEOMS, USER_OBSGEOMS,  & ! Inputs     !@@ 2p1
          N_USER_ANGLES, USER_ANGLES, N_USER_RELAZMS, USER_RELAZMS,       & ! Inputs
          FLUX_FACTOR, NBEAMS, BEAM_SZAS, EARTH_RADIUS, HEIGHT_GRID,      & ! Inputs
          DELTAU_INPUT, OMEGA_INPUT, ASYMM_INPUT, D2S_SCALING,            & ! Inputs
          THERMAL_BB_INPUT, LAMBERTIAN_ALBEDO, BRDF_F_0, BRDF_F, UBRDF_F, & ! Inputs
          EMISSIVITY, SURFBB, SLTERM_ISOTROPIC, SLTERM_F_0,               & ! Inputs  !@@ 2p3 (Sleave)
          INTENSITY_TOA, INTENSITY_BOA, FLUXES_TOA, FLUXES_BOA,           & ! Outputs !@@ 2p3 (Fluxes)
          RADLEVEL_UP, RADLEVEL_DN, N_GEOMETRIES,                         & ! Outputs
          STATUS_INPUTCHECK, C_NMESSAGES, C_MESSAGES, C_ACTIONS,          & ! Exception handling
          STATUS_EXECUTION,  E_MESSAGE, E_TRACE_1, E_TRACE_2 )              ! Exception handling

        call cpu_time(e2)
        exectime(thread) = exectime(thread) + e2-e1
        enddo

!  Exception handling

        IF ( STATUS_INPUTCHECK .eq. 1 ) THEN
          write(*,'(a,i4)')'INPUT Check failed from Baseline Run # ',T
          write(*,*)' - Number of Messages = ', C_NMESSAGES
          Do k = 1, C_NMESSAGES
            write(*,'(A,I3,A,A)')' - Message # ',K,': ', TRIM(C_MESSAGES(K))
            write(*,'(A,I3,A,A)')' - Action  # ',K,': ', TRIM(C_ACTIONS(K))
          ENDDO
          stop'Test_2S_only program aborted'
        ENDIF
        IF ( STATUS_EXECUTION .eq. 1 ) THEN
          write(*,'(a,i4)')'EXECUTION failed from Baseline Run # ',T
          write(*,*)' - Print 1 Message and 2 Traces'
          write(*,'(A)') TRIM(E_MESSAGE)
          write(*,'(A)') TRIM(E_TRACE_1)
          write(*,'(A)') TRIM(E_TRACE_2)
          stop'Test_2S_only program aborted'
        ENDIF

!  Save results

        INTENSITY_TOA_SAVE(1:n_geometries,T) = INTENSITY_TOA(1:n_geometries)
        INTENSITY_BOA_SAVE(1:n_geometries,T) = INTENSITY_BOA(1:n_geometries)
        FLUXES_TOA_SAVE(1:nbeams,1:2,t) = FLUXES_TOA(1:nbeams,1:2)
        FLUXES_BOA_SAVE(1:nbeams,1:2,t) = FLUXES_BOA(1:nbeams,1:2)

        if ( DO_2S_LEVELOUT ) Then
          do g = 1, n_geometries
            do n = 0, nlayers
               RADLEVEL_UP_SAVE(G,N,T) = RADLEVEL_UP(G,N)
               RADLEVEL_DN_SAVE(G,N,T) = RADLEVEL_DN(G,N)
            enddo
          enddo
        endif

!  End thread loop

      enddo

!  Timing tests
!  ============

! (skip if not set)

      if ( .not.do_timing ) go to 655

!  Open Existing file and read everything into execall,
!    Update first line, Add extra line with current timing checks

      k0 = n_user_obsgeoms
      if (      do_user_obsgeoms ) open(667,file='Timing_ObsGeom',status='unknown')
      if ( .not.do_user_obsgeoms ) open(667,file='Timing_Lattice',status='unknown')
      read(667,*)(nruns(n),n=1,6)
      nruns(k0)=nruns(k0)+1
      do n = 1, 6
         L = 0
         if ( n.eq.k0) then
            do m = 1, nruns(n) - 1
               L = L + 1 ; read(667,'(2i3,6f12.6)')ND,LD,(execall(L,n,t),t=1,nthreads)
            enddo
            l=l+1 ; execall(L,n,1:nthreads) = exectime(1:nthreads)
         else
            do m = 1, nruns(n)
               L = L + 1 ; read(667,'(2i3,6f12.6)')ND,LD,(execall(L,n,t),t=1,nthreads)
            enddo
         endif
      enddo
      close(667)

!  Open file again, and write again from fresh, including new and updated entries

      if (      do_user_obsgeoms ) open(667,file='Timing_ObsGeom',status='unknown')
      if ( .not.do_user_obsgeoms ) open(667,file='Timing_Lattice',status='unknown')
      write(667,'(6i2)')(nruns(n),n=1,6)
      do n = 1, 6
         do m = 1, nruns(n)
            write(667,'(2i3,6f12.6)')N,M,(execall(M,n,t),t=1,nthreads)
         enddo
      enddo
      close(667)

!  Screen output for this run

      if ( do_user_obsgeoms ) then
         do t = 1, nthreads
           write(*,'(a,i3,f12.6)')'User_Obsgeoms : thread # and time ',t,exectime(t)
         enddo
      else
         do t = 1, nthreads
           write(*,'(a,i3,f12.6)')'Lattice option: thread # and time ',t,exectime(t)
         enddo
      endif

!  Continuation point for skipping timing

655   continue

!  Prepare basic output file

      if ( do_user_obsgeoms)      OPEN(36,file = 'results_2s_tester.all_OBSGEOM', status = 'unknown')
      if ( .not.do_user_obsgeoms) OPEN(36,file = 'results_2s_tester.all_LATTICE', status = 'unknown')

!  Write results:

!  Intensity output

      write(36,'(/T32,a/T32,a/)')' 2S Intensity, threads 1-6', &
                                 '=========================='
      write(36,'(a,T31,6(a16,2x)/)')'Geometry    Level/Output', &
                                    '  Solar  No DM ','  Solar Yes DM ', &
                                    ' Thermal No DM ',' Thermal Yes DM', &
                                    'Solar+Thermal N','Solar+Thermal Y'
      do g = 1, n_geometries
        write(36,366)g,'Upwelling @ ',0.0d0, (INTENSITY_TOA_SAVE(G,T),T=1,NTHREADS)
        write(36,366)g,'Dnwelling @ ',23.0d0,(INTENSITY_BOA_SAVE(G,T),T=1,NTHREADS)
      enddo
366   format(i5,T11,a,f6.2,1x,6(1x,1pe16.8,1x))

!  Create output file (old format)

      do_old_outformat = .true.
      if (do_old_outformat) then

!  Mean output

      write(36,'(/T32,a/T32,a/)')'ACTINIC + REGULAR FLUXES, threads 1-2', &
                                 '====================================='
      write(36,'(a,T31,4(a16,2x)/)')' Sun SZA    Level/Output', &
            ' Actinic No DM ',' Actinic + DM  ',&
            ' Regular No DM ',' Regular + DM  '
      do v = 1, nbeams
         write(36,367)v,'Upwelling @ ',0.0d0, &
            (FLUXES_TOA_SAVE(v,1,t),t=1,2),(FLUXES_TOA_SAVE(v,2,t),t=1,2)
         write(36,367)v,'Dnwelling @ ',23.0d0, &
            (FLUXES_BOA_SAVE(v,1,t),t=1,2),(FLUXES_BOA_SAVE(v,2,t),t=1,2)
      enddo
367   format(i5,T11,a,f6.2,1x,4(1x,1pe16.8,1x))

      close(36)

!  Prepare level output file if selected

      if ( DO_2S_LEVELOUT ) Then
         if ( do_user_obsgeoms)      OPEN(36,file = 'results_2s_tester.LEVOUT_OBSGEOM', status = 'unknown')
         if ( .not.do_user_obsgeoms) OPEN(36,file = 'results_2s_tester.LEVOUT_LATTICE', status = 'unknown')

         write(36,'(/T32,a/T32,a/)')' 2S Intensity, threads 1-6', &
                                    '=========================='
         write(36,'(a,T31,6(a16,2x)/)')'Geometry    Level/Output', &
                                       '  Solar  No DM ','  Solar Yes DM ', &
                                       ' Thermal No DM ',' Thermal Yes DM', &
                                       'Solar+Thermal N','Solar+Thermal Y'
         do g = 1, n_geometries
            do n = 0, nlayers
               write(36,368)g,n,'Upwelling @ ',DBLE(N), (RADLEVEL_UP_SAVE(G,N,T),T=1,NTHREADS)
            enddo
            do n = 0, nlayers
               write(36,368)g,n,'Dnwelling @ ',DBLE(N), (RADLEVEL_DN_SAVE(G,N,T),T=1,NTHREADS)
            enddo
            write(36,*)
         enddo
368      format(i2,i3,T11,a,f6.2,2x,6(1x,1pe16.8,1x))

         close(36)
      endif

!  Create output file (new format)

      else

!  Mean output

      write(36,'(/T32,a/T32,a/)')'ACTINIC FLUXES, threads 1-6', &
                                 '==========================='
      write(36,'(a,T31,6(a16,2x)/)')' Sun SZA    Level/Output', &
                                    '  Solar  No DM ','  Solar Yes DM ', &
                                    ' Thermal No DM ',' Thermal Yes DM', &
                                    'Solar+Thermal N','Solar+Thermal Y'

      do v = 1, nbeams
         write(36,369)v,'Upwelling @ ',0.0d0 , (FLUXES_TOA_SAVE(v,1,t),t=1,NTHREADS)
         write(36,369)v,'Dnwelling @ ',23.0d0, (FLUXES_BOA_SAVE(v,1,t),t=1,NTHREADS)
      enddo

      write(36,'(/T32,a/T32,a/)')'REGULAR FLUXES, threads 1-6', &
                                 '==========================='
      write(36,'(a,T31,6(a16,2x)/)')' Sun SZA    Level/Output', &
                                    '  Solar  No DM ','  Solar Yes DM ', &
                                    ' Thermal No DM ',' Thermal Yes DM', &
                                    'Solar+Thermal N','Solar+Thermal Y'

      do v = 1, nbeams
         write(36,369)v,'Upwelling @ ',0.0d0 , (FLUXES_TOA_SAVE(v,2,t),t=1,NTHREADS)
         write(36,369)v,'Dnwelling @ ',23.0d0, (FLUXES_BOA_SAVE(v,2,t),t=1,NTHREADS)
      enddo
369   format(i5,T11,a,f6.2,1x,6(1x,1pe16.8,1x))

      close(36)

      if ( DO_2S_LEVELOUT ) Then
         if ( do_pentadiag_inverse ) then
            if ( do_user_obsgeoms)      OPEN(36,file = 'results_2s_tester_PDI.LEVOUT_OBSGEOM', status = 'unknown')
            if ( .not.do_user_obsgeoms) OPEN(36,file = 'results_2s_tester_PDI.LEVOUT_LATTICE', status = 'unknown')
         else
            if ( do_user_obsgeoms)      OPEN(36,file = 'results_2s_tester.LEVOUT_OBSGEOM', status = 'unknown')
            if ( .not.do_user_obsgeoms) OPEN(36,file = 'results_2s_tester.LEVOUT_LATTICE', status = 'unknown')
         endif

         write(36,'(/T32,a/T32,a/)')' 2S Intensity, threads 1-6', &
                                    '=========================='
         write(36,'(a,T31,6(a16,2x)/)')'Geometry    Level/Output', &
                                       '  Solar  No DM ','  Solar Yes DM ', &
                                       ' Thermal No DM ',' Thermal Yes DM', &
                                       'Solar+Thermal N','Solar+Thermal Y'
         do g = 1, n_geometries
            do n = 0, nlayers
               write(36,370)g,n,'Upwelling @',DBLE(N), (RADLEVEL_UP_SAVE(G,N,T),T=1,NTHREADS)
            enddo
            do n = 0, nlayers
               write(36,370)g,n,'Dnwelling @',DBLE(N), (RADLEVEL_DN_SAVE(G,N,T),T=1,NTHREADS)
            enddo
            write(36,*)
         enddo
370      format(i2,i3,T11,a,f6.2,2x,6(1x,1pe16.8,1x))

         close(36)
      endif

!  End format if block

      endif

!  Finish

      write(*,*)
      STOP 'successful run'

END PROGRAM test_2s_ONLY_2p4
