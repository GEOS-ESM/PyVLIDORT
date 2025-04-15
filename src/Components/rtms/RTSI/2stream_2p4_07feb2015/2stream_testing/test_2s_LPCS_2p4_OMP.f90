PROGRAM test_2s_LPCS_2p4_OMP

  USE twostream_ls_brdf_supplement_m
  USE twostream_lps_master_m
  USE twostream_lcs_master_m
  USE twostream_getPlanck

      implicit none

!  Notes 21 december 2012. Observational Geometry Inputs. Marked with !@@

!     Observation-Geometry New dimensioning.    MAX_USER_OBSGEOMS
!     Observation-Geometry input control.       DO_USER_OBSGEOMS
!     Observation-Geometry input control.       N_USER_OBSGEOMS
!     User-defined Observation Geometry angles. USER_OBSGEOMS

!  19 july    2013, ALL LEVEL output,      Version 2.2 (Tester 28 July)
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
      INTEGER, parameter :: one = 1.0d0

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

      INTEGER, PARAMETER :: max_atmoswfs   = 2
      INTEGER, PARAMETER :: max_surfacewfs = 2
      INTEGER, PARAMETER :: max_sleavewfs  = 1

      INTEGER, PARAMETER :: max_geometries = maxbeams * max_user_angles * max_user_relazms
      INTEGER, PARAMETER :: maxtotal       = 2*maxlayers

      INTEGER, PARAMETER :: max_omp_maxthreads = 2
      INTEGER, PARAMETER :: maxwvn = 10!1000

!  Directional Flags

      LOGICAL       :: DO_UPWELLING, DO_DNWELLING

!  Plane parallel and deltam-2stream scaling flags

      LOGICAL       :: DO_PLANE_PARALLEL, DO_D2S_SCALING

!  BRDF surface flag

      LOGICAL       :: DO_BRDF_SURFACE

!   !@@ Observational Geometry flag !@@

      LOGICAL       :: DO_USER_OBSGEOMS !@@

!  @@ Rob Spurr, 28 July 2013, Version 2.2, Levelout flag

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

!  Linearization flags. Add Sleave 2p3

      LOGICAL       :: DO_PROFILE_WFS, DO_COLUMN_WFS
      LOGICAL       :: DO_SURFACE_WFS, DO_SLEAVE_WFS

!  Linearization control. Add Sleave 2p3

      LOGICAL       :: LAYER_VARY_FLAG   ( MAXLAYERS )
      INTEGER       :: LAYER_VARY_NUMBER ( MAXLAYERS )
      INTEGER       :: N_COLUMN_WFS
      INTEGER       :: N_SURFACE_WFS, N_SLEAVE_WFS

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

      REAL(kind=dp) :: BEAM_SZAS    ( MAXBEAMS )
      REAL(kind=dp) :: USER_ANGLES  ( MAX_USER_ANGLES )
      REAL(kind=dp) :: USER_RELAZMS ( MAX_USER_RELAZMS )

!  Stream value

      REAL(kind=dp) :: STREAM_VALUE

!  Observational geometry input. [Same as LIDORT]. New 12/21/12 !@@

      INTEGER       :: N_USER_OBSGEOMS                    !@@
      REAL(kind=dp) :: USER_OBSGEOMS(MAX_USER_OBSGEOMS,3) !@@

!  Thermal inputs

      REAL(kind=dp) :: SURFBB
      REAL(kind=dp) :: THERMAL_BB_INPUT  ( 0:MAXLAYERS )

!  Lambertian Surface control (threaded)

      REAL(kind=dp) :: LAMBERTIAN_ALBEDO

!  BRDF Fourier components (NOT threaded)
!  0 and 1 Fourier components of BRDF, following order (same all threads)
!    incident solar directions,  reflected quadrature stream
!    incident quadrature stream, reflected quadrature stream
!    incident solar directions,  reflected user streams
!    incident quadrature stream, reflected user streams

      REAL(kind=dp) :: BRDF_F_0  ( 0:1, MAXBEAMS )
      REAL(kind=dp) :: BRDF_F    ( 0:1 )
!      REAL(kind=dp) :: UBRDF_F_0 ( 0:1, MAX_USER_ANGLES, MAXBEAMS )
      REAL(kind=dp) :: UBRDF_F   ( 0:1, MAX_USER_ANGLES )

!  Linearized BRDF fourier components
!  0 and 1 Fourier components of BRDF, following order (same all threads)
!    incident solar directions,  reflected quadrature stream
!    incident quadrature stream, reflected quadrature stream
!    incident solar directions,  reflected user streams
!    incident quadrature stream, reflected user streams

      REAL(kind=dp) :: LS_BRDF_F_0  ( MAX_SURFACEWFS, 0:1, MAXBEAMS )
      REAL(kind=dp) :: LS_BRDF_F    ( MAX_SURFACEWFS, 0:1 )
!     REAL(kind=dp) :: LS_UBRDF_F_0 ( MAX_SURFACEWFS, 0:1, MAX_USER_ANGLES, MAXBEAMS )
      REAL(kind=dp) :: LS_UBRDF_F   ( MAX_SURFACEWFS, 0:1, MAX_USER_ANGLES )

!  Emissivity

      REAL(kind=dp) :: EMISSIVITY
      REAL(kind=dp) :: LS_EMISSIVITY ( MAX_SURFACEWFS )

!  Version 2p3. 1/23/14. Introduce SLEAVE stuff
!    Do not require any first-order inputs (exact or Fourier)
!    Isotropic Surface leaving term (if flag set)
!  Fourier components of Surface-leaving terms:

      REAL(kind=dp) ::  SLTERM_ISOTROPIC ( MAXBEAMS )
      REAL(kind=dp) ::  SLTERM_F_0 ( 0:1, MAXBEAMS )

!  Exact Surface-Leaving term
!      REAL(kind=dp) ::  SLTERM_USERANGLES ( MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS )
!  Fourier components of Surface-leaving terms:
!    Every solar direction, SL-transmitted user streams. First order truncated
!      REAL(kind=dp) ::  USER_SLTERM_F_0 ( 0:1, MAX_USER_STREAMS, MAXBEAMS )
!  Input thread

!  Addition of SLEAVE WF inputs (Isotropic and Fourier components  diffuse-term)

      REAL(kind=dp)  :: LSSL_SLTERM_ISOTROPIC  ( MAX_SLEAVEWFS, MAXBEAMS )
      REAL(kind=dp)  :: LSSL_SLTERM_F_0        ( MAX_SLEAVEWFS, 0:1, MAXBEAMS )

!  Flux factor

      REAL(kind=dp) :: FLUX_FACTOR

!  height and earth radius

      REAL(kind=dp) :: EARTH_RADIUS
      REAL(kind=dp) :: HEIGHT_GRID ( 0:MAXLAYERS )

!  Atmospheric Optical properties

      REAL(kind=dp) :: DELTAU_INPUT( MAXLAYERS )
      REAL(kind=dp) :: OMEGA_INPUT ( MAXLAYERS )
      REAL(kind=dp) :: ASYMM_INPUT ( MAXLAYERS )
      REAL(kind=dp) :: D2S_SCALING ( MAXLAYERS )

!  Linearized optical properties

      REAL(kind=dp) :: L_DELTAU_INPUT( MAXLAYERS, MAX_ATMOSWFS )
      REAL(kind=dp) :: L_OMEGA_INPUT ( MAXLAYERS, MAX_ATMOSWFS )
      REAL(kind=dp) :: L_ASYMM_INPUT ( MAXLAYERS, MAX_ATMOSWFS )
      REAL(kind=dp) :: L_D2S_SCALING ( MAXLAYERS, MAX_ATMOSWFS )

!  Local Results (NOT THREADED NOW, Version 2.4)
!  -------------

!  Radiance & Radiance Jacobian output - TOA and BOA

      REAL(kind=dp) :: INTENSITY_TOA ( MAX_GEOMETRIES )
      REAL(kind=dp) :: INTENSITY_BOA ( MAX_GEOMETRIES )

      REAL(kind=dp) :: PROFILEWF_TOA ( MAX_GEOMETRIES, MAXLAYERS, MAX_ATMOSWFS )
      REAL(kind=dp) :: PROFILEWF_BOA ( MAX_GEOMETRIES, MAXLAYERS, MAX_ATMOSWFS )

      REAL(kind=dp) :: COLUMNWF_TOA  ( MAX_GEOMETRIES, MAX_ATMOSWFS )
      REAL(kind=dp) :: COLUMNWF_BOA  ( MAX_GEOMETRIES, MAX_ATMOSWFS )

      REAL(kind=dp) :: SURFACEWF_TOA ( MAX_GEOMETRIES, MAX_SURFACEWFS )
      REAL(kind=dp) :: SURFACEWF_BOA ( MAX_GEOMETRIES, MAX_SURFACEWFS )

      INTEGER       :: n_geometries

!  Radiance & Radiance Jacobian output - all levels
!     ! @@ Rob Spurr, 28 July 2013, Version 2.2 --> Optional Output at ALL LEVELS

      REAL(kind=dp) :: RADLEVEL_UP ( MAX_GEOMETRIES, 0:MAXLAYERS )
      REAL(kind=dp) :: RADLEVEL_DN ( MAX_GEOMETRIES, 0:MAXLAYERS )

      REAL(kind=dp) :: PROFJACLEVEL_UP ( MAX_GEOMETRIES, 0:MAXLAYERS, MAXLAYERS, MAX_ATMOSWFS )
      REAL(kind=dp) :: PROFJACLEVEL_DN ( MAX_GEOMETRIES, 0:MAXLAYERS, MAXLAYERS, MAX_ATMOSWFS )

      REAL(kind=dp) :: COLJACLEVEL_UP  ( MAX_GEOMETRIES, 0:MAXLAYERS, MAX_ATMOSWFS )
      REAL(kind=dp) :: COLJACLEVEL_DN  ( MAX_GEOMETRIES, 0:MAXLAYERS, MAX_ATMOSWFS )

      REAL(kind=dp) :: SURFJACLEVEL_UP ( MAX_GEOMETRIES, 0:MAXLAYERS, MAX_SURFACEWFS )
      REAL(kind=dp) :: SURFJACLEVEL_DN ( MAX_GEOMETRIES, 0:MAXLAYERS, MAX_SURFACEWFS )

!  Flux & Flux Jacobian output - TOA and BOA
!     ! @@ Rob Spurr, 05 November 2013, Version 2.3 --> Flux Output

      REAL(kind=dp) :: FLUXES_TOA ( MAXBEAMS, 2 )
      REAL(kind=dp) :: FLUXES_BOA ( MAXBEAMS, 2 )

      REAL(kind=dp) :: PROFJACFLUXES_TOA ( MAXBEAMS, 2, MAXLAYERS, MAX_ATMOSWFS )
      REAL(kind=dp) :: PROFJACFLUXES_BOA ( MAXBEAMS, 2, MAXLAYERS, MAX_ATMOSWFS )

      REAL(kind=dp) :: COLJACFLUXES_TOA  ( MAXBEAMS, 2, MAX_ATMOSWFS )
      REAL(kind=dp) :: COLJACFLUXES_BOA  ( MAXBEAMS, 2, MAX_ATMOSWFS )

      REAL(kind=dp) :: SURFJACFLUXES_TOA ( MAXBEAMS, 2, MAX_SURFACEWFS )
      REAL(kind=dp) :: SURFJACFLUXES_BOA ( MAXBEAMS, 2, MAX_SURFACEWFS )

!  Exception handling

!    1. Up to 100 Check Messages and actions

      INTEGER       :: STATUS_INPUTCHECK
      INTEGER       :: C_NMESSAGES
      CHARACTER*100 :: C_MESSAGES ( 0:MAXMESSAGES )
      CHARACTER*100 :: C_ACTIONS  ( 0:MAXMESSAGES )

!    2. Execution message and 2 Traces

      INTEGER       :: STATUS_EXECUTION
      CHARACTER*100 :: E_MESSAGE, E_TRACE_1, E_TRACE_2

!  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!          OTHER ARGUMENTS
!  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  External thread parameter

      INTEGER, PARAMETER :: maxthreads = 6


!  Saved Results, Baseline
!  -----------------------

!  Radiance & Radiance Jacobian output - TOA and BOA

      REAL(kind=dp) :: INTENSITY_TOA_BAS ( MAX_GEOMETRIES, MAXTHREADS, MAX_OMP_MAXTHREADS )
      REAL(kind=dp) :: INTENSITY_BOA_BAS ( MAX_GEOMETRIES, MAXTHREADS, MAX_OMP_MAXTHREADS )

      REAL(kind=dp) :: PROFILEWF_TOA_BAS ( MAX_GEOMETRIES, MAXLAYERS, MAX_ATMOSWFS, MAXTHREADS, MAX_OMP_MAXTHREADS )
      REAL(kind=dp) :: PROFILEWF_BOA_BAS ( MAX_GEOMETRIES, MAXLAYERS, MAX_ATMOSWFS, MAXTHREADS, MAX_OMP_MAXTHREADS )

      REAL(kind=dp) :: COLUMNWF_TOA_BAS ( MAX_GEOMETRIES, MAX_ATMOSWFS, MAXTHREADS, MAX_OMP_MAXTHREADS )
      REAL(kind=dp) :: COLUMNWF_BOA_BAS ( MAX_GEOMETRIES, MAX_ATMOSWFS, MAXTHREADS, MAX_OMP_MAXTHREADS )

      REAL(kind=dp) :: SURFACEWF_TOA_BAS ( MAX_GEOMETRIES, MAX_SURFACEWFS, MAXTHREADS, MAX_OMP_MAXTHREADS )
      REAL(kind=dp) :: SURFACEWF_BOA_BAS ( MAX_GEOMETRIES, MAX_SURFACEWFS, MAXTHREADS, MAX_OMP_MAXTHREADS )

      INTEGER       :: n_geometries_save

!  Radiance & Radiance Jacobian output - all levels
!     ! @@ Rob Spurr, 28 July 2013, Version 2.2 --> Optional Output at ALL LEVELS

      REAL(kind=dp) :: RADLEVEL_UP_BAS ( MAX_GEOMETRIES, 0:MAXLAYERS, MAXTHREADS, MAX_OMP_MAXTHREADS )
      REAL(kind=dp) :: RADLEVEL_DN_BAS ( MAX_GEOMETRIES, 0:MAXLAYERS, MAXTHREADS, MAX_OMP_MAXTHREADS )

      REAL(kind=dp) :: PROFJACLEVEL_UP_BAS ( MAX_GEOMETRIES, 0:MAXLAYERS, MAXLAYERS, MAX_ATMOSWFS, MAXTHREADS, MAX_OMP_MAXTHREADS )
      REAL(kind=dp) :: PROFJACLEVEL_DN_BAS ( MAX_GEOMETRIES, 0:MAXLAYERS, MAXLAYERS, MAX_ATMOSWFS, MAXTHREADS, MAX_OMP_MAXTHREADS )

      REAL(kind=dp) :: COLJACLEVEL_UP_BAS ( MAX_GEOMETRIES, 0:MAXLAYERS, MAX_ATMOSWFS, MAXTHREADS, MAX_OMP_MAXTHREADS )
      REAL(kind=dp) :: COLJACLEVEL_DN_BAS ( MAX_GEOMETRIES, 0:MAXLAYERS, MAX_ATMOSWFS, MAXTHREADS, MAX_OMP_MAXTHREADS )

      REAL(kind=dp) :: SURFJACLEVEL_UP_BAS ( MAX_GEOMETRIES, 0:MAXLAYERS, MAX_SURFACEWFS, MAXTHREADS, MAX_OMP_MAXTHREADS )
      REAL(kind=dp) :: SURFJACLEVEL_DN_BAS ( MAX_GEOMETRIES, 0:MAXLAYERS, MAX_SURFACEWFS, MAXTHREADS, MAX_OMP_MAXTHREADS )

!  Flux & Flux Jacobian output - TOA and BOA
!     ! @@ Rob Spurr, 05 November 2013, Version 2.3 --> Flux Output

      REAL(kind=dp) :: FLUXES_TOA_BAS ( MAXBEAMS, 2, MAXTHREADS, MAX_OMP_MAXTHREADS )
      REAL(kind=dp) :: FLUXES_BOA_BAS ( MAXBEAMS, 2, MAXTHREADS, MAX_OMP_MAXTHREADS )

      REAL(kind=dp) :: PROFJACFLUXES_TOA_BAS ( MAXBEAMS, 2, MAXLAYERS, MAX_ATMOSWFS, MAXTHREADS, MAX_OMP_MAXTHREADS )
      REAL(kind=dp) :: PROFJACFLUXES_BOA_BAS ( MAXBEAMS, 2, MAXLAYERS, MAX_ATMOSWFS, MAXTHREADS, MAX_OMP_MAXTHREADS )

      REAL(kind=dp) :: COLJACFLUXES_TOA_BAS ( MAXBEAMS, 2, MAX_ATMOSWFS, MAXTHREADS, MAX_OMP_MAXTHREADS )
      REAL(kind=dp) :: COLJACFLUXES_BOA_BAS ( MAXBEAMS, 2, MAX_ATMOSWFS, MAXTHREADS, MAX_OMP_MAXTHREADS )

      REAL(kind=dp) :: SURFJACFLUXES_TOA_BAS ( MAXBEAMS, 2, MAX_SURFACEWFS, MAXTHREADS, MAX_OMP_MAXTHREADS )
      REAL(kind=dp) :: SURFJACFLUXES_BOA_BAS ( MAXBEAMS, 2, MAX_SURFACEWFS, MAXTHREADS, MAX_OMP_MAXTHREADS )


!  Saved Results, Perturbed values for FD Jacobians
!  ------------------------------------------------

!  Radiance & Radiance Jacobian output - TOA and BOA

      REAL(kind=dp) :: INTENSITY_TOA_PROF_PT ( MAX_GEOMETRIES, MAXLAYERS, MAX_ATMOSWFS, MAXTHREADS, MAX_OMP_MAXTHREADS )
      REAL(kind=dp) :: INTENSITY_BOA_PROF_PT ( MAX_GEOMETRIES, MAXLAYERS, MAX_ATMOSWFS, MAXTHREADS, MAX_OMP_MAXTHREADS )

      REAL(kind=dp) :: INTENSITY_TOA_COL_PT ( MAX_GEOMETRIES, MAX_ATMOSWFS, MAXTHREADS, MAX_OMP_MAXTHREADS )
      REAL(kind=dp) :: INTENSITY_BOA_COL_PT ( MAX_GEOMETRIES, MAX_ATMOSWFS, MAXTHREADS, MAX_OMP_MAXTHREADS )

      REAL(kind=dp) :: INTENSITY_TOA_SURF_PT ( MAX_GEOMETRIES, MAX_SURFACEWFS, MAXTHREADS, MAX_OMP_MAXTHREADS )
      REAL(kind=dp) :: INTENSITY_BOA_SURF_PT ( MAX_GEOMETRIES, MAX_SURFACEWFS, MAXTHREADS, MAX_OMP_MAXTHREADS )

!  Radiance & Radiance Jacobian output - all levels
!     ! @@ Rob Spurr, 28 July 2013, Version 2.2 --> Optional Output at ALL LEVELS

      REAL(kind=dp) :: RADLEVEL_UP_PROF_PT ( MAX_GEOMETRIES, 0:MAXLAYERS, MAXLAYERS, MAX_ATMOSWFS, MAXTHREADS, MAX_OMP_MAXTHREADS )
      REAL(kind=dp) :: RADLEVEL_DN_PROF_PT ( MAX_GEOMETRIES, 0:MAXLAYERS, MAXLAYERS, MAX_ATMOSWFS, MAXTHREADS, MAX_OMP_MAXTHREADS )

      REAL(kind=dp) :: RADLEVEL_UP_COL_PT ( MAX_GEOMETRIES, 0:MAXLAYERS, MAX_ATMOSWFS, MAXTHREADS, MAX_OMP_MAXTHREADS )
      REAL(kind=dp) :: RADLEVEL_DN_COL_PT ( MAX_GEOMETRIES, 0:MAXLAYERS, MAX_ATMOSWFS, MAXTHREADS, MAX_OMP_MAXTHREADS )

      REAL(kind=dp) :: RADLEVEL_UP_SURF_PT ( MAX_GEOMETRIES, 0:MAXLAYERS, MAX_SURFACEWFS, MAXTHREADS, MAX_OMP_MAXTHREADS )
      REAL(kind=dp) :: RADLEVEL_DN_SURF_PT ( MAX_GEOMETRIES, 0:MAXLAYERS, MAX_SURFACEWFS, MAXTHREADS, MAX_OMP_MAXTHREADS )

!  Flux & Flux Jacobian output - TOA and BOA
!     ! @@ Rob Spurr, 05 November 2013, Version 2.3 --> Flux Output

      REAL(kind=dp) :: FLUXES_TOA_PROF_PT ( MAXBEAMS, 2, MAXLAYERS, MAX_ATMOSWFS, MAXTHREADS, MAX_OMP_MAXTHREADS )
      REAL(kind=dp) :: FLUXES_BOA_PROF_PT ( MAXBEAMS, 2, MAXLAYERS, MAX_ATMOSWFS, MAXTHREADS, MAX_OMP_MAXTHREADS )

      REAL(kind=dp) :: FLUXES_TOA_COL_PT ( MAXBEAMS, 2, MAX_ATMOSWFS, MAXTHREADS, MAX_OMP_MAXTHREADS )
      REAL(kind=dp) :: FLUXES_BOA_COL_PT ( MAXBEAMS, 2, MAX_ATMOSWFS, MAXTHREADS, MAX_OMP_MAXTHREADS )

      REAL(kind=dp) :: FLUXES_TOA_SURF_PT ( MAXBEAMS, 2, MAX_SURFACEWFS, MAXTHREADS, MAX_OMP_MAXTHREADS )
      REAL(kind=dp) :: FLUXES_BOA_SURF_PT ( MAXBEAMS, 2, MAX_SURFACEWFS, MAXTHREADS, MAX_OMP_MAXTHREADS )


!  Other flags

      LOGICAL       :: DO_FULLQUADRATURE

!  Help variables

      INTEGER       :: f,n,n6,ndum,ldum,t,k,i,q,v
      INTEGER       :: nthreads, thread
      REAL(kind=dp) :: kd, gaer, waer, taer, parcel, raywt, aerwt
      REAL(kind=dp) :: aersca, aerext, molsca, totsca, totext, c0, molabs
      REAL(kind=dp) :: molomg(maxlayers),molext(maxlayers), m1, m2
      REAL(kind=dp) :: raymoms ( 0:2, maxlayers ), omega, column
      REAL(kind=dp) :: eps, epsfac, brdf_par
      REAL(kind=dp) :: albedo_save
      REAL(kind=dp) :: windspeed_save

!  Thermal emission setups

      LOGICAL           :: THERMFAIL
      CHARACTER(LEN=70) :: THERMALMESSAGE
      REAL(kind=dp)     :: SURFTEMP, AIRTEMPS ( 0:MAXLAYERS )
      REAL(kind=dp)     :: WNUMLO, WNUMHI
      INTEGER           :: NTEMPS, SMALLV

!  Cox-Munk Surface control (threaded)

      REAL(kind=dp) :: WINDSPEED ( MAXTHREADS )

!  BRDF inputs

      INTEGER       :: NSTREAMS_BRDF
      LOGICAL       :: LAMBERTIAN_KERNEL_FLAG ( MAX_BRDF_KERNELS )
      LOGICAL       :: DO_SHADOW_EFFECT
      INTEGER       :: N_BRDF_KERNELS
      INTEGER       :: WHICH_BRDF   ( MAX_BRDF_KERNELS )
      REAL(kind=dp) :: BRDF_FACTORS ( MAX_BRDF_KERNELS )
      INTEGER       :: N_BRDF_PARAMETERS ( MAX_BRDF_KERNELS )
      REAL(kind=dp) :: BRDF_PARAMETERS   ( MAX_BRDF_KERNELS, MAX_BRDF_PARAMETERS )

      LOGICAL       :: DO_KERNEL_FACTOR_WFS ( MAX_BRDF_KERNELS )
      LOGICAL       :: DO_KERNEL_PARAMS_WFS ( MAX_BRDF_KERNELS, MAX_BRDF_PARAMETERS )

      INTEGER       :: N_KERNEL_FACTOR_WFS
      INTEGER       :: N_KERNEL_PARAMS_WFS
      LOGICAL       :: DO_KPARAMS_DERIVS ( MAX_BRDF_PARAMETERS )

!  BRDF Exception handling. !@@ Added, 12/31/12

      INTEGER       :: status_brdfsup
      CHARACTER*100 :: message_brdf, action_brdf

!  OpenMP tests (general)

      INTEGER       :: TEST, N_OMP_TESTS, OMP_MAXTHREADS
      INTEGER       :: TID, OMP_NTHREADS
      INTEGER       :: wvn, numwvn, ompt

      CHARACTER (LEN=1)  :: charTID
      CHARACTER (LEN=1)  :: CH1 
      CHARACTER (LEN=5)  :: CH2,CH3
      CHARACTER (LEN=30) :: TAIL

!  OpenMP tests (timing)

      INTEGER       :: n_core, time_divider
      REAL          :: omp_e1, omp_e2

!  OpenMP functions

      INTEGER :: OMP_GET_NUM_THREADS, &
                 OMP_GET_THREAD_NUM

!  Set test parameters

      N_CORE = 2
      N_OMP_TESTS = 2

      numwvn   = maxwvn
      nthreads = 6

!  Begin test loop

      DO TEST = 1, N_OMP_TESTS
      !DO TEST = 2, N_OMP_TESTS
      !DO TEST = 1, 1
      !DO TEST = 2, 2

      write(*,*)
      write(*,*) '******************************************'
      write(*,'(1x,a,i1)') 'Doing test ',TEST

!  Set total number of threads (optional)
!      Note: if not set, OpenMP will set the number of threads
!            to the number of cores on the system
!      Note: if not set, OpenMP will default to dividing up the
!            wvn loop below equally among the available threads

      if ( test == 1 ) then
        OMP_MAXTHREADS = 1
      elseif ( test == 2 ) then
        OMP_MAXTHREADS = 2
      endif

      !write(*,*)
      !write(*,*) 'OMP_MAXTHREADS = ',OMP_MAXTHREADS

!  Set total number of threads

      !write(*,*)
      !write(*,*) 'Setting total number of threads'
      CALL OMP_SET_NUM_THREADS(OMP_MAXTHREADS)

!  Start timing

      call cpu_time(omp_e1)

!  Set some control logicals

      DO_UPWELLING = .true.
      DO_DNWELLING = .true.

!  Control read

      open(44,file = '2s_2p4_driver.inp',status = 'old')
      read(44,*) DO_SOLAR_SOURCES
      read(44,*) DO_THERMAL_EMISSION
      read(44,*) DO_SURFACE_EMISSION
      read(44,*) DO_BRDF_SURFACE
      read(44,*) DO_FULLQUADRATURE
      read(44,*) DO_D2S_SCALING
      read(44,*) DO_PLANE_PARALLEL
      read(44,*) DO_USER_OBSGEOMS     ! @@ 2p1
      read(44,*)
      read(44,*)
      read(44,*) DO_2S_LEVELOUT       ! @@ 2p2
      read(44,*) DO_MVOUT_ONLY        ! @@ 2p3
      read(44,*) DO_ADDITIONAL_MVOUT  ! @@ 2p3
      read(44,*) DO_PENTADIAG_INVERSE ! @@ 2p3
      close(44)

!  Modify test file input

!      DO_PLANE_PARALLEL = .false.
      DO_BRDF_SURFACE   = .true.!.false.

!  Define other control inputs not in the input file

      DO_SURFACE_LEAVING = .false.
      DO_SL_ISOTROPIC    = .false.
      DO_SLEAVE_WFS      = .false.
      N_SLEAVE_WFS       = 0

!  Taylor control

      TAYLOR_SMALL = 1.0d-03
      TAYLOR_ORDER = 3

!  Scale factor and BVP Index are pre-set

      BVPINDEX       = 1
      BVPSCALEFACTOR = 1.0d0

!  Control integers

      n_user_obsgeoms = max_user_obsgeoms
      nbeams          = n_user_obsgeoms
      n_user_angles   = n_user_obsgeoms
      n_user_relazms  = n_user_obsgeoms

! mick fix 1/16/2015 - save the max value of N_GEOMETRIES for the current setup
!                      (note: not necessarily equal to the max parameter value
!                             set at the top of the driver)
      if (.not. DO_USER_OBSGEOMS) then
        n_geometries_save = nbeams * n_user_angles * n_user_relazms
      else
        n_geometries_save = n_user_obsgeoms
      end if

      nlayers  = 23
      ntotal   = 2 * nlayers

      nstreams_brdf  = 50

!  Geophysical factors

      flux_factor  = 1.0d0
      earth_radius = 6371.0d0

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

!  Get the pre-prepared atmosphere

      open(45,file='data/input_atmos.dat',status='old' )
      read(45,'(i5,1p25e18.9)')ldum, (raymoms(0,n),n=1,nlayers)
      read(45,'(i5,1p25e18.9)')ldum, (raymoms(1,n),n=1,nlayers)
      !Note: raymoms(2,n) has already been multiplied by 5 = 2*L+1 where L=2
      read(45,'(i5,1p25e18.9)')ldum, (raymoms(2,n),n=1,nlayers)
      height_grid(0) = 60.0d0
      do n = 1, nlayers
         read(45,'(i4,f12.5,1p6e16.7)')ndum,height_grid(n), &
            molext(n),molomg(n),kd,kd,kd,kd
      enddo
      close(45)
      height_grid(nlayers+1:) = 0.0d0

!  Total column absorption

      c0 = 0.0d0
      do n = 1, nlayers
        molabs = molext(n) * ( 1.0d0 - molomg(n) )
        c0 = c0 + molabs
      enddo

!  Add Aerosols bottom 6 layers, spread evenly

      n6 = nlayers - 6; gaer = 0.8d0 ; waer = 0.95d0 ; taer = 0.5d0

!  Thermal stuff
!  -------------

!  Array of temperatures should be 24 levels

      open(1,file= 'data/input_temp23.dat', status='old')
      read(1,*)SURFTEMP
      read(1,*)NTEMPS
      do n = 0, NTEMPS
         read(1,*)ndum, AIRTEMPS(NTEMPS-N)
      enddo
      close(1)

!  Surface stuff

      albedo_save    =  0.2d0
      windspeed_save = 10.0d0

!  Two choices of stream value................ CHOOSE One !!!!!

      if ( do_fullquadrature ) then
        stream_value = dsqrt(1.0d0 / 3.0d0 )
      else
        stream_value = 0.5d0
      endif

!  Begin file open loop

      do ompt=1,OMP_MAXTHREADS

!  Open output files

        write(ch1,'(i1)')   omp_maxthreads
        write(ch2,'(i5.5)') numwvn 
        write(ch3,'(i5.5)') ompt*numwvn/omp_maxthreads
        tail = '_nt' // ch1 // '_nwn' // ch2 // '_wvn' // ch3

      if ( do_user_obsgeoms ) then
         if ( do_pentadiag_inverse ) then
           OPEN(10*ompt+36,file = 'results_2s_LPCS_tester_PDI.TOAUP_OBSGEOM' // trim(tail), status = 'unknown')
           OPEN(10*ompt+37,file = 'results_2s_LPCS_tester_PDI.BOADN_OBSGEOM' // trim(tail), status = 'unknown')
         else
           OPEN(10*ompt+36,file = 'results_2s_LPCS_tester.TOAUP_OBSGEOM' // trim(tail), status = 'unknown')
           OPEN(10*ompt+37,file = 'results_2s_LPCS_tester.BOADN_OBSGEOM' // trim(tail), status = 'unknown')
         endif
      else
         if ( do_pentadiag_inverse ) then
           OPEN(10*ompt+36,file = 'results_2s_LPCS_tester_PDI.TOAUP_LATTICE' // trim(tail), status = 'unknown')
           OPEN(10*ompt+37,file = 'results_2s_LPCS_tester_PDI.BOADN_LATTICE' // trim(tail), status = 'unknown')
         else
           OPEN(10*ompt+36,file = 'results_2s_LPCS_tester.TOAUP_LATTICE' // trim(tail), status = 'unknown')
           OPEN(10*ompt+37,file = 'results_2s_LPCS_tester.BOADN_LATTICE' // trim(tail), status = 'unknown')
         endif
      endif

      if ( do_2S_LEVELOUT ) then
         if ( do_user_obsgeoms ) then
            if ( do_pentadiag_inverse ) then
               OPEN(10*ompt+336,file = 'results_2s_LPCS_tester_PDI.LEVOUT_UP_OBSGEOM' // trim(tail), status = 'unknown')
               OPEN(10*ompt+337,file = 'results_2s_LPCS_tester_PDI.LEVOUT_DN_OBSGEOM' // trim(tail), status = 'unknown')
            else
               OPEN(10*ompt+336,file = 'results_2s_LPCS_tester.LEVOUT_UP_OBSGEOM' // trim(tail), status = 'unknown')
               OPEN(10*ompt+337,file = 'results_2s_LPCS_tester.LEVOUT_DN_OBSGEOM' // trim(tail), status = 'unknown')
            endif
         else
            if ( do_pentadiag_inverse ) then
               OPEN(10*ompt+336,file = 'results_2s_LPCS_tester_PDI.LEVOUT_UP_LATTICE' // trim(tail), status = 'unknown')
               OPEN(10*ompt+337,file = 'results_2s_LPCS_tester_PDI.LEVOUT_DN_LATTICE' // trim(tail), status = 'unknown')
            else
               OPEN(10*ompt+336,file = 'results_2s_LPCS_tester.LEVOUT_UP_LATTICE' // trim(tail), status = 'unknown')
               OPEN(10*ompt+337,file = 'results_2s_LPCS_tester.LEVOUT_DN_LATTICE' // trim(tail), status = 'unknown')
            endif
         endif
      endif

!  End file open loop

      enddo

!  Define epsilon for FD perturbations 

      eps = 1.0d-03

!  Only master thread does this

      !OMP_NTHREADS = OMP_GET_NUM_THREADS()
      !write(*,*)
      !write(*,'(1x,a,i1)') 'Total number of threads (to start) = ', OMP_NTHREADS

!      go to 656

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!                   Part 1: PROFILE/SURFACE WFS TESTS                 @
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  Initialize saved output for Intensity and associated Jacobians

      INTENSITY_TOA_BAS     = 0.0d0
      INTENSITY_TOA_PROF_PT = 0.0d0
      INTENSITY_TOA_SURF_PT = 0.0d0
      PROFILEWF_TOA_BAS     = 0.0d0
      SURFACEWF_TOA_BAS     = 0.0d0

      INTENSITY_BOA_BAS     = 0.0d0
      INTENSITY_BOA_PROF_PT = 0.0d0
      INTENSITY_BOA_SURF_PT = 0.0d0
      PROFILEWF_BOA_BAS     = 0.0d0
      SURFACEWF_BOA_BAS     = 0.0d0

!  Initialize saved output for Level Radiances and associated Jacobians

      RADLEVEL_UP_BAS       = 0.0d0
      RADLEVEL_UP_PROF_PT   = 0.0d0
      RADLEVEL_UP_SURF_PT   = 0.0d0
      PROFJACLEVEL_UP_BAS   = 0.0d0
      SURFJACLEVEL_UP_BAS   = 0.0d0

      RADLEVEL_DN_BAS       = 0.0d0
      RADLEVEL_DN_PROF_PT   = 0.0d0
      RADLEVEL_DN_SURF_PT   = 0.0d0
      PROFJACLEVEL_DN_BAS   = 0.0d0
      SURFJACLEVEL_DN_BAS   = 0.0d0

!  Initialize saved output for Fluxes and associated Jacobians

      FLUXES_TOA_BAS        = 0.0d0
      FLUXES_TOA_PROF_PT    = 0.0d0
      FLUXES_TOA_SURF_PT    = 0.0d0
      PROFJACFLUXES_TOA_BAS = 0.0d0
      SURFJACFLUXES_TOA_BAS = 0.0d0

      FLUXES_BOA_BAS        = 0.0d0
      FLUXES_BOA_PROF_PT    = 0.0d0
      FLUXES_BOA_SURF_PT    = 0.0d0
      PROFJACFLUXES_BOA_BAS = 0.0d0
      SURFJACFLUXES_BOA_BAS = 0.0d0

!     Begin parallel region

!$OMP PARALLEL IF(OMP_MAXTHREADS > 1) &
!$OMP   DEFAULT (PRIVATE) &
!$OMP   SHARED (OMP_NTHREADS, &
!$OMP     DO_UPWELLING, DO_DNWELLING, &
!$OMP     DO_SURFACE_LEAVING, DO_SL_ISOTROPIC, DO_SLEAVE_WFS, N_SLEAVE_WFS, &
!$OMP     BVPINDEX, BVPSCALEFACTOR, TAYLOR_ORDER, TAYLOR_SMALL, &
!$OMP     NUMWVN, NTHREADS, NLAYERS, NTOTAL, &
!$OMP     N_USER_OBSGEOMS, N_USER_ANGLES, NSTREAMS_BRDF, &
!$OMP     BEAM_SZAS, USER_ANGLES, USER_RELAZMS, USER_OBSGEOMS, &
!$OMP     STREAM_VALUE, FLUX_FACTOR, EARTH_RADIUS, &
!$OMP     HEIGHT_GRID, MOLEXT, MOLOMG, RAYMOMS, N6, GAER, WAER, TAER, &
!$OMP     NTEMPS, AIRTEMPS, SURFTEMP, &
!$OMP     WINDSPEED_SAVE, ALBEDO_SAVE, BRDF_PAR, eps, &
!$OMP     INTENSITY_TOA_BAS, INTENSITY_BOA_BAS, &
!$OMP     PROFILEWF_TOA_BAS, PROFILEWF_BOA_BAS, INTENSITY_TOA_PROF_PT, INTENSITY_BOA_PROF_PT, &
!$OMP     SURFACEWF_TOA_BAS, SURFACEWF_BOA_BAS, INTENSITY_TOA_SURF_PT, INTENSITY_BOA_SURF_PT, &
!$OMP     RADLEVEL_UP_BAS, RADLEVEL_DN_BAS, & 
!$OMP     PROFJACLEVEL_UP_BAS, PROFJACLEVEL_DN_BAS, RADLEVEL_UP_PROF_PT, RADLEVEL_DN_PROF_PT, &
!$OMP     SURFJACLEVEL_UP_BAS, SURFJACLEVEL_DN_BAS, RADLEVEL_UP_SURF_PT, RADLEVEL_DN_SURF_PT, &
!$OMP     FLUXES_TOA_BAS, FLUXES_BOA_BAS, &
!$OMP     PROFJACFLUXES_TOA_BAS, PROFJACFLUXES_BOA_BAS, FLUXES_TOA_PROF_PT, FLUXES_BOA_PROF_PT, &
!$OMP     SURFJACFLUXES_TOA_BAS, SURFJACFLUXES_BOA_BAS, FLUXES_TOA_SURF_PT, FLUXES_BOA_SURF_PT, &
!$OMP     N_GEOMETRIES_SAVE) &
!$OMP   FIRSTPRIVATE (DO_SOLAR_SOURCES, DO_THERMAL_EMISSION, DO_SURFACE_EMISSION, &
!$OMP     DO_BRDF_SURFACE, DO_FULLQUADRATURE, DO_D2S_SCALING, &
!$OMP     DO_PLANE_PARALLEL, DO_USER_OBSGEOMS, DO_2S_LEVELOUT, &
!$OMP     DO_MVOUT_ONLY, DO_ADDITIONAL_MVOUT, DO_PENTADIAG_INVERSE, &
!$OMP     NBEAMS, N_USER_RELAZMS)

!     Obtain thread number

      TID = OMP_GET_THREAD_NUM()
      !TID = 0
      
!     Obtain and display total number of threads and local thread

      IF (TID == 0) THEN
        OMP_NTHREADS = OMP_GET_NUM_THREADS()
        !OMP_NTHREADS = 1
        !write(*,*)
        !write(*,'(1x,a,i1)') 'Total number of threads (inside parallel region) = ', OMP_NTHREADS
        write(*,'(1x,a,i1,a)') 'Running 2S driver with ',OMP_NTHREADS,' thread(s)'
        write(*,*)
      END IF

      write(*,'(1x,a,i1,a)') 'Thread ',TID,' beginning its wavenumber loop'

!$OMP DO

!  Begin "wavenumber" loop

      do wvn = 1, numwvn

      !write(*,'(1x,a,i1,a,i5)') 'Thread  = ',TID,' wvn = ',wvn
      write(*,*)
      write(*,*)'Doing OMP thread ',TID,' Profile WF baseline: Intensities, 2 Profile WFs + 1 Surface WF'
      write(*,*) 'wvn = ',wvn

!  Initialize some inputs

      deltau_input = 0.0d0
      omega_input  = 0.0d0
      asymm_input  = 0.0d0
      d2s_scaling  = 0.0d0

      BRDF_F_0 = 0.0d0
      BRDF_F   = 0.0d0
      UBRDF_F  = 0.0d0

      lambertian_albedo = 0.0d0

      L_deltau_input = 0.0d0
      L_omega_input  = 0.0d0
      L_asymm_input  = 0.0d0
      L_d2s_scaling  = 0.0d0

      LS_BRDF_F_0 = 0.0d0
      LS_BRDF_F   = 0.0d0
      LS_UBRDF_F  = 0.0d0

!  Control for 2 profile WFs

      do_profile_wfs = .true.
      do_column_wfs  = .false.
      do_surface_wfs = .true.
      do n = 1, nlayers
         layer_vary_flag(n)   = .true.
         layer_vary_number(n) = 2
      enddo
      n_column_wfs   = 0
      n_surface_wfs  = 1

!  Control for SLEAVE

      DO_SURFACE_LEAVING = .FALSE.
      DO_SL_ISOTROPIC    = .FALSE.
      SLTERM_ISOTROPIC   = 0.0D0
      SLTERM_F_0         = 0.0D0
      DO_SLEAVE_WFS      = .FALSE.
      N_SLEAVE_WFS       = 0
      LSSL_SLTERM_ISOTROPIC = 0.0D0
      LSSL_SLTERM_F_0       = 0.0D0

!  Baseline calculation 1 : RADIANCE + PROFILE/SURFACE WFS
!  =======================================================

!  Start thread loop

      do thread = 1, nthreads
!      do thread = 3, 3
        t = thread

!  Thread 1: 2S, Solar only   , no   delta-M scaling
!  Thread 2: 2S, Solar only   , with delta-M scaling
!  Thread 3: 2S, Thermal only , no   delta-M scaling
!  Thread 4: 2S, Thermal only , with delta-M scaling
!  Thread 5: 2S, Solar+Thermal, no   delta-M scaling
!  Thread 6: 2S, Solar+Thermal, with delta-M scaling

        call thread_setter &
          ( thread, n_user_obsgeoms, & 
            do_solar_sources, do_thermal_emission, &
            do_d2s_scaling, do_surface_emission,   &
            nbeams, n_user_relazms )
        !write(*,*)'Doing thread # ', thread
        !write(*,'(1x,a,i1,a,i1)') 'Doing OMP thread ',TID,' LIDORT thread ',thread

        column = 1.0d0

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
          !emissivity    =  1.0d0 - albedo_save
          !ls_emissivity = -1.0d0
        else
          surfbb        = 0.0d0
          emissivity    = 0.0d0
          ls_emissivity = 0.0d0
        end if

!  Create optical properties, every thread

        do n = 1, n6
          deltau_input(n) = molext(n)
          omega_input(n)  = molomg(n)
          asymm_input(n)  = 0.0d0

          molsca = molomg(n) * molext(n)
          molabs = molext(n) - molsca
          m1     = raymoms(2,n) * molsca ; m2 = 0.0d0
          d2s_scaling(n)  = ( m1 + m2 ) / molsca / 5.0d0


          !1st linearized oprop set (partially normalized) -
          !  par=1 wrt mol abs tau
          omega = omega_input(n)
          l_deltau_input(n,1) = molabs
          l_omega_input (n,1) = - molabs * omega/deltau_input(n)
        enddo

        parcel = taer / ( height_grid(n6) - height_grid(nlayers) )
        do n = n6 + 1, nlayers
          aerext = Parcel * ( height_grid(n-1) - height_grid(n) )
          molsca = molomg(n) * molext(n)
          molabs = molext(n) - molsca

!          molabs = molabs * 1.001d0
!          aerext = aerext * 1.001d0

          aersca = aerext * waer
          totext = molsca + molabs + aerext
          totsca = molsca + aersca
          raywt  = molsca / totsca
          aerwt  = aersca / totsca
          omega  = totsca / totext

          deltau_input(n) = totext
          omega_input(n)  = omega
          asymm_input(n)  = gaer * aerwt

          m1 = raymoms(2,n) * molsca
          m2 = 5.0d0 * gaer * gaer * aersca
          d2s_scaling(n)  = ( m1 + m2 ) / totsca / 5.0d0


          !1st linearized oprop set (partially normalized) -
          !  par=1 wrt mol abs tau
          l_deltau_input(n,1) = molabs
          l_omega_input (n,1) = - molabs * omega/deltau_input(n)


          !2nd linearized oprop set (partially normalized) -
          !  par=2 wrt aero ext tau

! @@@@@@@@@@@@@@@@@@@@@@
!  This set was wrong
!          l_deltau_input(n,2) = aerext
!          l_omega_input(n,2)  = aerext * (1.0d0 - omega) / deltau_input(n)
!          l_asymm_input(n,2)  = aerext * (gaer - asymm_input(n)) / totsca
!          l_d2s_scaling(n,2)  = aerext * (gaer*gaer - d2s_scaling(n)) / totsca

!  Rob fixed, 14 nov 12
          l_deltau_input(n,2) = aerext
          l_omega_input(n,2)  = aerext * (waer - omega) / deltau_input(n)
          l_asymm_input(n,2)  = aerwt  * (gaer - asymm_input(n))
          l_d2s_scaling(n,2)  = aerwt  * (gaer*gaer - d2s_scaling(n))
! @@@@@@@@@@@@@@@@@@@@@@

!        if ( n.gt.18)write(*,*)n,l_asymm_input(n,2,1),asymm_input(n,1)
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

          !Linearized BRDF supplement inputs
          DO_KERNEL_FACTOR_WFS      = .FALSE.
          DO_KERNEL_PARAMS_WFS      = .FALSE.

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
          !      TWOSTREAM_LS_BRDFMASTER using the DO_SHADOW_EFFECT flag

          !Linearized BRDF supplement inputs
          DO_KERNEL_PARAMS_WFS(1,1) = .TRUE.

          !Save BRDF slope parameter for perturbation comparison below
          BRDF_PAR = BRDF_PARAMETERS(1,1)

          !Call linearized BRDF supplement
          CALL TWOSTREAM_LS_BRDFMASTER &
            ( MAXBEAMS, MAX_USER_ANGLES, MAX_USER_OBSGEOMS,      & ! Dimensions !@@
              MAXSTREAMS_BRDF, MAX_BRDF_KERNELS,                 & ! Dimensions
              MAX_BRDF_PARAMETERS, MAX_SURFACEWFS,               & ! Dimensions
              DO_SOLAR_SOURCES, DO_USER_OBSGEOMS,                & ! Inputs !@@
              LAMBERTIAN_KERNEL_FLAG,                            & ! Inputs
              DO_SHADOW_EFFECT, DO_SURFACE_EMISSION,             & ! Inputs
              NBEAMS, N_USER_ANGLES, N_USER_OBSGEOMS,            & ! Inputs !@@
              BEAM_SZAS, USER_ANGLES, USER_OBSGEOMS,             & ! Inputs !@@
              STREAM_VALUE, NSTREAMS_BRDF,                       & ! Inputs
              N_BRDF_KERNELS, WHICH_BRDF, BRDF_FACTORS,          & ! Inputs
              N_BRDF_PARAMETERS, BRDF_PARAMETERS,                & ! Inputs
              DO_KERNEL_FACTOR_WFS, DO_KERNEL_PARAMS_WFS,        & ! Inputs
              DO_KPARAMS_DERIVS, N_SURFACE_WFS,                  & ! Outputs
              N_KERNEL_FACTOR_WFS, N_KERNEL_PARAMS_WFS,          & ! Outputs
              BRDF_F_0, BRDF_F, UBRDF_F, EMISSIVITY,             & ! Outputs
              LS_BRDF_F_0, LS_BRDF_F, LS_UBRDF_F, LS_EMISSIVITY, & ! Outputs
              STATUS_BRDFSUP, MESSAGE_BRDF, ACTION_BRDF )          ! Outputs

        else
          !Lambertian surface

          BRDF_PAR = albedo_save
          LAMBERTIAN_ALBEDO = albedo_save
          emissivity    =  1.0d0 - albedo_save
          ls_emissivity = -1.0d0
        end if

!  Exception handling on BRDF

        if ( DO_BRDF_SURFACE ) then
           if ( STATUS_BRDFSUP .ne. 0 ) then
              write(*,'(a,i4)')'BRDF Lin supplement Check failed from Baseline Profile WF Run # ',T
              write(*,*)' - Print 1 Message and 1 Action'
              write(*,'(A)') TRIM(MESSAGE_BRDF)
              write(*,'(A)') TRIM(ACTION_BRDF)
              stop'Test_2S_LPCS program aborted'
           endif
        endif

!  Call to Linearized model

        CALL TWOSTREAM_LPS_MASTER &
        ( MAXLAYERS, MAXTOTAL, MAXMESSAGES, MAXBEAMS, MAX_GEOMETRIES,      & ! Dimensions
          MAX_USER_ANGLES, MAX_USER_RELAZMS, MAX_USER_OBSGEOMS,            & ! Dimensions !@@ 2p1
          MAX_ATMOSWFS, MAX_SURFACEWFS, MAX_SLEAVEWFS,                     & ! Dimensions !@@ 2p3 (Add Sleave)
          DO_UPWELLING, DO_DNWELLING, DO_PLANE_PARALLEL, DO_2S_LEVELOUT,   & ! Inputs     !@@ 2p2
          DO_MVOUT_ONLY, DO_ADDITIONAL_MVOUT,                              & ! Inputs     !@@ 2p3
          DO_SOLAR_SOURCES, DO_THERMAL_EMISSION, DO_SURFACE_EMISSION,      & ! Inputs
          DO_D2S_SCALING, DO_BRDF_SURFACE, DO_USER_OBSGEOMS,               & ! Inputs     !@@ 2p1
          DO_SURFACE_LEAVING, DO_SL_ISOTROPIC, DO_PENTADIAG_INVERSE,       & ! Input !@@ 2p3   6/25/14
          BVPINDEX, BVPSCALEFACTOR, TAYLOR_ORDER, TAYLOR_SMALL,            & ! Input !@@ 2p3/4 6/25/14, 1/7/15
          NLAYERS, NTOTAL, STREAM_VALUE, N_USER_OBSGEOMS, USER_OBSGEOMS,   & ! Inputs     !@@ 2p1
          N_USER_ANGLES, USER_ANGLES, N_USER_RELAZMS, USER_RELAZMS,        & ! Inputs
          FLUX_FACTOR, NBEAMS, BEAM_SZAS, EARTH_RADIUS, HEIGHT_GRID,       & ! Inputs
          DELTAU_INPUT, OMEGA_INPUT, ASYMM_INPUT, D2S_SCALING,             & ! Inputs
          THERMAL_BB_INPUT, LAMBERTIAN_ALBEDO, BRDF_F_0, BRDF_F, UBRDF_F,  & ! Inputs
          EMISSIVITY, SURFBB, SLTERM_ISOTROPIC, SLTERM_F_0,                & ! Inputs  !@@ 2p3 (Add Sleave)
          DO_PROFILE_WFS, DO_SURFACE_WFS, DO_SLEAVE_WFS,                   & ! Inputs  !@@ 2p3 (Add Sleave)
          LAYER_VARY_FLAG, LAYER_VARY_NUMBER, N_SURFACE_WFS, N_SLEAVE_WFS, & ! Inputs
          LSSL_SLTERM_ISOTROPIC, LSSL_SLTERM_F_0,                          & ! Inputs  !@@ 2p3 (New sleave)
          L_DELTAU_INPUT, L_OMEGA_INPUT, L_ASYMM_INPUT, L_D2S_SCALING,     & ! Inputs
          LS_BRDF_F_0, LS_BRDF_F, LS_UBRDF_F, LS_EMISSIVITY,               & ! Inputs
          INTENSITY_TOA, PROFILEWF_TOA, SURFACEWF_TOA,                     & ! In/Out
          INTENSITY_BOA, PROFILEWF_BOA, SURFACEWF_BOA,                     & ! In/Out
          RADLEVEL_UP, RADLEVEL_DN,  N_GEOMETRIES,                         & ! Outputs !@@ 2p2
          PROFJACLEVEL_UP,PROFJACLEVEL_DN,SURFJACLEVEL_UP,SURFJACLEVEL_DN, & ! Outputs !@@ 2p2
          FLUXES_TOA, PROFJACFLUXES_TOA, SURFJACFLUXES_TOA,                & ! Outputs !@@ 2p3
          FLUXES_BOA, PROFJACFLUXES_BOA, SURFJACFLUXES_BOA,                & ! Outputs !@@ 2p3
          STATUS_INPUTCHECK, C_NMESSAGES, C_MESSAGES, C_ACTIONS,           & ! Exception handling
          STATUS_EXECUTION,  E_MESSAGE, E_TRACE_1, E_TRACE_2 )               ! Exception handling

!  Exception handling

        IF ( STATUS_INPUTCHECK .eq. 1 ) THEN
          write(*,'(a,i4)')'INPUT Check failed from Baseline Profile WF Run # ',T
          write(*,*)' - Number of Messages = ', C_NMESSAGES
          DO k = 1, C_NMESSAGES
            write(*,'(A,I3,A,A)')' - Message # ',K,': ', TRIM(C_MESSAGES(K))
            write(*,'(A,I3,A,A)')' - Action  # ',K,': ', TRIM(C_ACTIONS(K))
          ENDDO
          stop 'Test_2S_LPCS program aborted'
        ENDIF
        IF ( STATUS_EXECUTION .eq. 1 ) THEN
          write(*,'(a,i4)')'EXECUTION failed from Baseline Profile WF Run # ',T
          write(*,*)' - Print 1 Message and 2 Traces'
          write(*,'(A)') TRIM(E_MESSAGE)
          write(*,'(A)') TRIM(E_TRACE_1)
          write(*,'(A)') TRIM(E_TRACE_2)
          stop 'Test_2S_LPCS program aborted'
        ENDIF

!  Save the last results produced by the OMP thread

        IF (wvn == (TID+1)*numwvn/OMP_NTHREADS) THEN

!if (t == 1) write(*,*) 'wvn inside = ',wvn

          INTENSITY_TOA_BAS ( 1:N_GEOMETRIES, T, TID+1 ) = INTENSITY_TOA ( 1:N_GEOMETRIES )
          INTENSITY_BOA_BAS ( 1:N_GEOMETRIES, T, TID+1 ) = INTENSITY_BOA ( 1:N_GEOMETRIES )

          PROFILEWF_TOA_BAS ( 1:N_GEOMETRIES, 1:NLAYERS, :, T, TID+1 ) = &
              PROFILEWF_TOA ( 1:N_GEOMETRIES, 1:NLAYERS, : )
          PROFILEWF_BOA_BAS ( 1:N_GEOMETRIES, 1:NLAYERS, :, T, TID+1 ) = &
              PROFILEWF_BOA ( 1:N_GEOMETRIES, 1:NLAYERS, : )

          SURFACEWF_TOA_BAS ( 1:N_GEOMETRIES, 1:N_SURFACE_WFS, T, TID+1 ) = &
              SURFACEWF_TOA ( 1:N_GEOMETRIES, 1:N_SURFACE_WFS )
          SURFACEWF_BOA_BAS ( 1:N_GEOMETRIES, 1:N_SURFACE_WFS, T, TID+1 ) = &
              SURFACEWF_BOA ( 1:N_GEOMETRIES, 1:N_SURFACE_WFS )

          IF ( do_2s_LEVELOUT ) THEN
            RADLEVEL_UP_BAS ( 1:N_GEOMETRIES, 0:NLAYERS, T, TID+1 ) = &
                RADLEVEL_UP ( 1:N_GEOMETRIES, 0:NLAYERS )
            RADLEVEL_DN_BAS ( 1:N_GEOMETRIES, 0:NLAYERS, T, TID+1 ) = &
                RADLEVEL_DN ( 1:N_GEOMETRIES, 0:NLAYERS )

            PROFJACLEVEL_UP_BAS ( 1:N_GEOMETRIES, 0:NLAYERS, 1:NLAYERS, :, T, TID+1 ) = &
                PROFJACLEVEL_UP ( 1:N_GEOMETRIES, 0:NLAYERS, 1:NLAYERS, : )
            PROFJACLEVEL_DN_BAS ( 1:N_GEOMETRIES, 0:NLAYERS, 1:NLAYERS, :, T, TID+1 ) = &
                PROFJACLEVEL_DN ( 1:N_GEOMETRIES, 0:NLAYERS, 1:NLAYERS, : )

            SURFJACLEVEL_UP_BAS ( 1:N_GEOMETRIES, 0:NLAYERS, 1:N_SURFACE_WFS, T, TID+1 ) = &
                SURFJACLEVEL_UP ( 1:N_GEOMETRIES, 0:NLAYERS, 1:N_SURFACE_WFS )
            SURFJACLEVEL_DN_BAS ( 1:N_GEOMETRIES, 0:NLAYERS, 1:N_SURFACE_WFS, T, TID+1 ) = &
                SURFJACLEVEL_DN ( 1:N_GEOMETRIES, 0:NLAYERS, 1:N_SURFACE_WFS )
          ENDIF

          FLUXES_TOA_BAS ( 1:NBEAMS, :, T, TID+1 ) = &
              FLUXES_TOA ( 1:NBEAMS, : )
          FLUXES_BOA_BAS ( 1:NBEAMS, :, T, TID+1 ) = &
              FLUXES_BOA ( 1:NBEAMS, : )

          PROFJACFLUXES_TOA_BAS ( 1:NBEAMS, :, 1:NLAYERS, :, T, TID+1 ) = &
              PROFJACFLUXES_TOA ( 1:NBEAMS, :, 1:NLAYERS, : )
          PROFJACFLUXES_BOA_BAS ( 1:NBEAMS, :, 1:NLAYERS, :, T, TID+1 ) = &
              PROFJACFLUXES_BOA ( 1:NBEAMS, :, 1:NLAYERS, : )

          SURFJACFLUXES_TOA_BAS ( 1:NBEAMS, :, 1:N_SURFACE_WFS, T, TID+1 ) = &
              SURFJACFLUXES_TOA ( 1:NBEAMS, :, 1:N_SURFACE_WFS )
          SURFJACFLUXES_BOA_BAS ( 1:NBEAMS, :, 1:N_SURFACE_WFS, T, TID+1 ) = &
              SURFJACFLUXES_BOA ( 1:NBEAMS, :, 1:N_SURFACE_WFS )
        END IF

!  End thread loop

      enddo

!  Preparing for FD Profile Jacobian calculations
!  ==============================================

!  Zeroing

      L_deltau_input = 0.0d0
      L_omega_input  = 0.0d0
      L_asymm_input  = 0.0d0
      L_d2s_scaling  = 0.0d0

      LS_BRDF_F_0    = 0.0d0
      LS_BRDF_F      = 0.0d0
      LS_UBRDF_F     = 0.0d0
      LS_EMISSIVITY  = 0.0d0

!  Control for no WFs

      do_profile_wfs = .false.
      do_surface_wfs = .false.
      do n = 1, nlayers
         layer_vary_flag(n)   = .false.
         layer_vary_number(n) = 0
      enddo
      n_column_wfs   = 0
      n_surface_wfs  = 0

!  FD perturbation

      !eps = 1.0d-03
      epsfac = 1.0d0 + eps


!  FD calculation 1 : Perturbation on Cox-Munk slope parameter
!  ===========================================================

      !write(*,*)
      !if (do_brdf_surface) then
      !  write(*,*) 'Doing OMP thread ',TID,'.  Doing FD calculation 1: perturb Cox-Munk slope parameter'
      !else
      !  write(*,*) 'Doing OMP thread ',TID,'.  Doing FD calculation 1: perturb Lambertian albedo'
      !endif

!  Start thread loop

      do thread = 1, nthreads
!      do thread = 3, 3
        t = thread

!  Thread 1: 2S, Solar only   , no   delta-M scaling
!  Thread 2: 2S, Solar only   , with delta-M scaling
!  Thread 3: 2S, Thermal only , no   delta-M scaling
!  Thread 4: 2S, Thermal only , with delta-M scaling
!  Thread 5: 2S, Solar+Thermal, no   delta-M scaling
!  Thread 6: 2S, Solar+Thermal, with delta-M scaling

        call thread_setter &
          ( thread, n_user_obsgeoms, & 
            do_solar_sources, do_thermal_emission, &
            do_d2s_scaling, do_surface_emission,   &
            nbeams, n_user_relazms )

        !write(*,*)'Doing thread # ', thread
        !write(*,'(1x,a,i1,a,i1)') 'Doing OMP thread ',TID,' LIDORT thread ',thread

        column = 1.0d0

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
          !emissivity    =  1.0d0 - albedo_save
        else
          surfbb        = 0.0d0
          emissivity    = 0.0d0
        end if

!  Create optical properties, every thread

        do n = 1, n6
          deltau_input(n) = molext(n)
          omega_input(n)  = molomg(n)
          asymm_input(n)  = 0.0d0

          molsca = molomg(n) * molext(n)
          molabs = molext(n) - molsca
          m1     = raymoms(2,n) * molsca ; m2 = 0.0d0
          d2s_scaling(n)  = ( m1 + m2 ) / molsca / 5.0d0
        enddo

        parcel = taer / ( height_grid(n6) - height_grid(nlayers) )
        do n = n6 + 1, nlayers
          aerext = Parcel * ( height_grid(n-1) - height_grid(n) )
          molsca = molomg(n) * molext(n)
          molabs = molext(n) - molsca

!          molabs = molabs * 1.001d0
!          aerext = aerext * 1.001d0

          aersca = aerext * waer
          totext = molsca + molabs + aerext
          totsca = molsca + aersca
          raywt  = molsca / totsca
          aerwt  = aersca / totsca
          omega  = totsca / totext

          deltau_input(n) = totext
          omega_input(n)  = omega
          asymm_input(n)  = gaer * aerwt

          m1 = raymoms(2,n) * molsca
          m2 = 5.0d0 * gaer * gaer * aersca
          d2s_scaling(n)  = ( m1 + m2 ) / totsca / 5.0d0
        enddo

!  Set up surface BRDF

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

          !Linearized BRDF supplement inputs
          DO_KERNEL_FACTOR_WFS      = .FALSE.
          DO_KERNEL_PARAMS_WFS      = .FALSE.

          !(2) Now, for the Cox-Munk BRDF, define:

          !Basic BRDF supplement inputs
          N_BRDF_KERNELS            = 1

          WHICH_BRDF(1)             = 9 !Cox-Munk
          BRDF_FACTORS(1)           = 1.0D0
          N_BRDF_PARAMETERS(1)      = 3
          WINDSPEED(T)              = windspeed_save
          BRDF_PARAMETERS(1,1)      = (0.003d0+0.00512d0*WINDSPEED(T))*epsfac !Slope squared
          BRDF_PARAMETERS(1,2)      = 1.33D0*1.33D0 !Square of refractive index
          !Note: The Cox-Munk shadow effect parameter (the 3rd
          !      parameter BRDF_PARAMETERS(1,3)) is set internally by
          !      TWOSTREAM_LS_BRDFMASTER using the DO_SHADOW_EFFECT flag

          !Call linearized BRDF supplement
          CALL TWOSTREAM_LS_BRDFMASTER &
            ( MAXBEAMS, MAX_USER_ANGLES, MAX_USER_OBSGEOMS,      & ! Dimensions !@@
              MAXSTREAMS_BRDF, MAX_BRDF_KERNELS,                 & ! Dimensions
              MAX_BRDF_PARAMETERS, MAX_SURFACEWFS,               & ! Dimensions
              DO_SOLAR_SOURCES, DO_USER_OBSGEOMS,                & ! Inputs !@@
              LAMBERTIAN_KERNEL_FLAG,                            & ! Inputs
              DO_SHADOW_EFFECT, DO_SURFACE_EMISSION,             & ! Inputs
              NBEAMS, N_USER_ANGLES, N_USER_OBSGEOMS,            & ! Inputs !@@
              BEAM_SZAS, USER_ANGLES, USER_OBSGEOMS,             & ! Inputs !@@
              STREAM_VALUE, NSTREAMS_BRDF,                       & ! Inputs
              N_BRDF_KERNELS, WHICH_BRDF, BRDF_FACTORS,          & ! Inputs
              N_BRDF_PARAMETERS, BRDF_PARAMETERS,                & ! Inputs
              DO_KERNEL_FACTOR_WFS, DO_KERNEL_PARAMS_WFS,        & ! Inputs
              DO_KPARAMS_DERIVS, N_SURFACE_WFS,                  & ! Outputs
              N_KERNEL_FACTOR_WFS, N_KERNEL_PARAMS_WFS,          & ! Outputs
              BRDF_F_0, BRDF_F, UBRDF_F, EMISSIVITY,             & ! Outputs
              LS_BRDF_F_0, LS_BRDF_F, LS_UBRDF_F, LS_EMISSIVITY, & ! Outputs
              STATUS_BRDFSUP, MESSAGE_BRDF, ACTION_BRDF )          ! Outputs
        else
          !Lambertian surface

          LAMBERTIAN_ALBEDO = albedo_save * epsfac
          emissivity = 1.0d0 - LAMBERTIAN_ALBEDO
        end if

!  Exception handling on BRDF

        if ( DO_BRDF_SURFACE ) then
           if ( STATUS_BRDFSUP .ne. 0 ) then
              write(*,'(a,i4)')'BRDF Lin supplement Check failed from Perturbed Profile WF Run # ',T
              write(*,*)' - Print 1 Message and 1 Action'
              write(*,'(A)') TRIM(MESSAGE_BRDF)
              write(*,'(A)') TRIM(ACTION_BRDF)
              stop'Test_2S_LPCS program aborted'
           endif
        endif

!  Call to Linearized model

        CALL TWOSTREAM_LPS_MASTER &
        ( MAXLAYERS, MAXTOTAL, MAXMESSAGES, MAXBEAMS, MAX_GEOMETRIES,      & ! Dimensions
          MAX_USER_ANGLES, MAX_USER_RELAZMS, MAX_USER_OBSGEOMS,            & ! Dimensions !@@ 2p1
          MAX_ATMOSWFS, MAX_SURFACEWFS, MAX_SLEAVEWFS,                     & ! Dimensions !@@ 2p3 (Add Sleave)
          DO_UPWELLING, DO_DNWELLING, DO_PLANE_PARALLEL, DO_2S_LEVELOUT,   & ! Inputs     !@@ 2p2
          DO_MVOUT_ONLY, DO_ADDITIONAL_MVOUT,                              & ! Inputs     !@@ 2p3
          DO_SOLAR_SOURCES, DO_THERMAL_EMISSION, DO_SURFACE_EMISSION,      & ! Inputs
          DO_D2S_SCALING, DO_BRDF_SURFACE, DO_USER_OBSGEOMS,               & ! Inputs     !@@ 2p1
          DO_SURFACE_LEAVING, DO_SL_ISOTROPIC, DO_PENTADIAG_INVERSE,       & ! Input !@@ 2p3   6/25/14
          BVPINDEX, BVPSCALEFACTOR, TAYLOR_ORDER, TAYLOR_SMALL,            & ! Input !@@ 2p3/4 6/25/14, 1/7/15
          NLAYERS, NTOTAL, STREAM_VALUE, N_USER_OBSGEOMS, USER_OBSGEOMS,   & ! Inputs     !@@ 2p1
          N_USER_ANGLES, USER_ANGLES, N_USER_RELAZMS, USER_RELAZMS,        & ! Inputs
          FLUX_FACTOR, NBEAMS, BEAM_SZAS, EARTH_RADIUS, HEIGHT_GRID,       & ! Inputs
          DELTAU_INPUT, OMEGA_INPUT, ASYMM_INPUT, D2S_SCALING,             & ! Inputs
          THERMAL_BB_INPUT, LAMBERTIAN_ALBEDO, BRDF_F_0, BRDF_F, UBRDF_F,  & ! Inputs
          EMISSIVITY, SURFBB, SLTERM_ISOTROPIC, SLTERM_F_0,                & ! Inputs  !@@ 2p3 (Add Sleave)
          DO_PROFILE_WFS, DO_SURFACE_WFS, DO_SLEAVE_WFS,                   & ! Inputs  !@@ 2p3 (Add Sleave)
          LAYER_VARY_FLAG, LAYER_VARY_NUMBER, N_SURFACE_WFS, N_SLEAVE_WFS, & ! Inputs
          LSSL_SLTERM_ISOTROPIC, LSSL_SLTERM_F_0,                          & ! Inputs  !@@ 2p3 (New sleave)
          L_DELTAU_INPUT, L_OMEGA_INPUT, L_ASYMM_INPUT, L_D2S_SCALING,     & ! Inputs
          LS_BRDF_F_0, LS_BRDF_F, LS_UBRDF_F, LS_EMISSIVITY,               & ! Inputs
          INTENSITY_TOA, PROFILEWF_TOA, SURFACEWF_TOA,                     & ! In/Out
          INTENSITY_BOA, PROFILEWF_BOA, SURFACEWF_BOA,                     & ! In/Out
          RADLEVEL_UP, RADLEVEL_DN,  N_GEOMETRIES,                         & ! Outputs !@@ 2p2
          PROFJACLEVEL_UP,PROFJACLEVEL_DN,SURFJACLEVEL_UP,SURFJACLEVEL_DN, & ! Outputs !@@ 2p2
          FLUXES_TOA, PROFJACFLUXES_TOA, SURFJACFLUXES_TOA,                & ! Outputs !@@ 2p3
          FLUXES_BOA, PROFJACFLUXES_BOA, SURFJACFLUXES_BOA,                & ! Outputs !@@ 2p3
          STATUS_INPUTCHECK, C_NMESSAGES, C_MESSAGES, C_ACTIONS,           & ! Exception handling
          STATUS_EXECUTION,  E_MESSAGE, E_TRACE_1, E_TRACE_2 )               ! Exception handling

!  Exception handling

        IF ( STATUS_INPUTCHECK .eq. 1 ) THEN
          write(*,'(a,i4)')'INPUT Check failed from Perturbed Profile WF Run # ',T
          write(*,*)' - Number of Messages = ', C_NMESSAGES
          DO k = 1, C_NMESSAGES
            write(*,'(A,I3,A,A)')' - Message # ',K,': ', TRIM(C_MESSAGES(K))
            write(*,'(A,I3,A,A)')' - Action  # ',K,': ', TRIM(C_ACTIONS(K))
          ENDDO
          stop'Test_2S_LPCS program aborted'
        ENDIF
        IF ( STATUS_EXECUTION .eq. 1 ) THEN
          write(*,'(a,i4)')'EXECUTION failed from Perturbed Profile WF Run # ',T
          write(*,*)' - Print 1 Message and 2 Traces'
          write(*,'(A)') TRIM(E_MESSAGE)
          write(*,'(A)') TRIM(E_TRACE_1)
          write(*,'(A)') TRIM(E_TRACE_2)
          stop'Test_2S_LPCS program aborted'
        ENDIF

!  Save the last results produced by the OMP thread

        IF (wvn == (TID+1)*numwvn/OMP_NTHREADS) THEN
          Q = 1
          INTENSITY_TOA_SURF_PT(1:N_GEOMETRIES,Q,T,TID+1) = INTENSITY_TOA(1:N_GEOMETRIES)
          INTENSITY_BOA_SURF_PT(1:N_GEOMETRIES,Q,T,TID+1) = INTENSITY_BOA(1:N_GEOMETRIES)

          if ( do_2s_LEVELOUT ) THEN
            RADLEVEL_UP_SURF_PT(1:N_GEOMETRIES,0:NLAYERS,Q,T,TID+1) = RADLEVEL_UP(1:N_GEOMETRIES,0:NLAYERS)
            RADLEVEL_DN_SURF_PT(1:N_GEOMETRIES,0:NLAYERS,Q,T,TID+1) = RADLEVEL_DN(1:N_GEOMETRIES,0:NLAYERS)
          endif

          FLUXES_TOA_SURF_PT(1:NBEAMS,:,Q,T,TID+1) = FLUXES_TOA(1:NBEAMS,:)
          FLUXES_BOA_SURF_PT(1:NBEAMS,:,Q,T,TID+1) = FLUXES_BOA(1:NBEAMS,:)
        END IF

!  End thread loop

      enddo


!  FD calculation 2 : Perturbation on mol abs tau
!  ==============================================

      !write(*,*)
      !write(*,*)'Doing OMP thread ',TID,' Doing FD calculation 2: perturb mol abs tau'

!  Start thread loop

      do thread = 1, nthreads
!      do thread = 3, 3
        t = thread

!  Thread 1: 2S, Solar only   , no   delta-M scaling
!  Thread 2: 2S, Solar only   , with delta-M scaling
!  Thread 3: 2S, Thermal only , no   delta-M scaling
!  Thread 4: 2S, Thermal only , with delta-M scaling
!  Thread 5: 2S, Solar+Thermal, no   delta-M scaling
!  Thread 6: 2S, Solar+Thermal, with delta-M scaling

        call thread_setter &
          ( thread, n_user_obsgeoms, & 
            do_solar_sources, do_thermal_emission, &
            do_d2s_scaling, do_surface_emission,   &
            nbeams, n_user_relazms )
        !write(*,*)'Doing thread # ', thread
        !write(*,'(1x,a,i1,a,i1)') 'Doing OMP thread ',TID,' LIDORT thread ',thread

        column = 1.0d0

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
          !emissivity    =  1.0d0 - albedo_save
        else
          surfbb        = 0.0d0
          emissivity    = 0.0d0
        end if

!  Start layer loop

        do i=1,nlayers

!  Create optical properties, every thread

          do n = 1, n6
            molsca = molomg(n) * molext(n)
            molabs = molext(n) - molsca

            if (n == i) then
              deltau_input(n) = molsca + molabs*epsfac
            else
              deltau_input(n) = molsca + molabs
            endif
            omega_input(n)  = molsca/deltau_input(n)
            asymm_input(n)  = 0.0d0

            m1 = raymoms(2,n) * molsca ; m2 = 0.0d0
            d2s_scaling(n)  = ( m1 + m2 ) / molsca / 5.0d0
          enddo

          parcel = taer / ( height_grid(n6) - height_grid(nlayers) )
          do n = n6 + 1, nlayers
            molsca = molomg(n) * molext(n)
            molabs = molext(n) - molsca

            aerext = Parcel * ( height_grid(n-1) - height_grid(n) )
            aersca = aerext * waer

            if (n == i) then
              totext = molsca + molabs*epsfac + aerext
            else
              totext = molsca + molabs + aerext
            endif
            totsca = molsca + aersca
            raywt  = molsca / totsca
            aerwt  = aersca / totsca
            omega  = totsca / totext

            deltau_input(n) = totext
            omega_input(n)  = omega
            asymm_input(n)  = gaer * aerwt

            m1 = raymoms(2,n) * molsca
            m2 = 5.0d0 * gaer * gaer * aersca
            d2s_scaling(n)  = ( m1 + m2 ) / totsca / 5.0d0
          enddo

!  Set up surface BRDF

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

            !Linearized BRDF supplement inputs
            DO_KERNEL_FACTOR_WFS      = .FALSE.
            DO_KERNEL_PARAMS_WFS      = .FALSE.

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
            !      TWOSTREAM_LS_BRDFMASTER using the DO_SHADOW_EFFECT flag

            !Call linearized BRDF supplement
            CALL TWOSTREAM_LS_BRDFMASTER &
            ( MAXBEAMS, MAX_USER_ANGLES, MAX_USER_OBSGEOMS,      & ! Dimensions !@@
              MAXSTREAMS_BRDF, MAX_BRDF_KERNELS,                 & ! Dimensions
              MAX_BRDF_PARAMETERS, MAX_SURFACEWFS,               & ! Dimensions
              DO_SOLAR_SOURCES, DO_USER_OBSGEOMS,                & ! Inputs !@@
              LAMBERTIAN_KERNEL_FLAG,                            & ! Inputs
              DO_SHADOW_EFFECT, DO_SURFACE_EMISSION,             & ! Inputs
              NBEAMS, N_USER_ANGLES, N_USER_OBSGEOMS,            & ! Inputs !@@
              BEAM_SZAS, USER_ANGLES, USER_OBSGEOMS,             & ! Inputs !@@
              STREAM_VALUE, NSTREAMS_BRDF,                       & ! Inputs
              N_BRDF_KERNELS, WHICH_BRDF, BRDF_FACTORS,          & ! Inputs
              N_BRDF_PARAMETERS, BRDF_PARAMETERS,                & ! Inputs
              DO_KERNEL_FACTOR_WFS, DO_KERNEL_PARAMS_WFS,        & ! Inputs
              DO_KPARAMS_DERIVS, N_SURFACE_WFS,                  & ! Outputs
              N_KERNEL_FACTOR_WFS, N_KERNEL_PARAMS_WFS,          & ! Outputs
              BRDF_F_0, BRDF_F, UBRDF_F, EMISSIVITY,             & ! Outputs
              LS_BRDF_F_0, LS_BRDF_F, LS_UBRDF_F, LS_EMISSIVITY, & ! Outputs
              STATUS_BRDFSUP, MESSAGE_BRDF, ACTION_BRDF )          ! Outputs
          else
            !Lambertian surface

            LAMBERTIAN_ALBEDO = albedo_save
            emissivity =  1.0d0 - albedo_save
          end if

!  Exception handling on BRDF

          if ( DO_BRDF_SURFACE ) then
             if ( STATUS_BRDFSUP .ne. 0 ) then
                write(*,'(a,i4)')'BRDF Lin supplement Check failed from Perturbed Profile WF Run # ',T
                write(*,*)' - Print 1 Message and 1 Action'
                write(*,'(A)') TRIM(MESSAGE_BRDF)
                write(*,'(A)') TRIM(ACTION_BRDF)
                stop'Test_2S_LPCS program aborted'
             endif
          endif

!  Call to Linearized model

          CALL TWOSTREAM_LPS_MASTER &
          ( MAXLAYERS, MAXTOTAL, MAXMESSAGES, MAXBEAMS, MAX_GEOMETRIES,      & ! Dimensions
            MAX_USER_ANGLES, MAX_USER_RELAZMS, MAX_USER_OBSGEOMS,            & ! Dimensions !@@ 2p1
            MAX_ATMOSWFS, MAX_SURFACEWFS, MAX_SLEAVEWFS,                     & ! Dimensions !@@ 2p3 (Add Sleave)
            DO_UPWELLING, DO_DNWELLING, DO_PLANE_PARALLEL, DO_2S_LEVELOUT,   & ! Inputs     !@@ 2p2
            DO_MVOUT_ONLY, DO_ADDITIONAL_MVOUT,                              & ! Inputs     !@@ 2p3
            DO_SOLAR_SOURCES, DO_THERMAL_EMISSION, DO_SURFACE_EMISSION,      & ! Inputs
            DO_D2S_SCALING, DO_BRDF_SURFACE, DO_USER_OBSGEOMS,               & ! Inputs     !@@ 2p1
            DO_SURFACE_LEAVING, DO_SL_ISOTROPIC, DO_PENTADIAG_INVERSE,       & ! Input !@@ 2p3   6/25/14
            BVPINDEX, BVPSCALEFACTOR, TAYLOR_ORDER, TAYLOR_SMALL,            & ! Input !@@ 2p3/4 6/25/14, 1/7/15
            NLAYERS, NTOTAL, STREAM_VALUE, N_USER_OBSGEOMS, USER_OBSGEOMS,   & ! Inputs     !@@ 2p1
            N_USER_ANGLES, USER_ANGLES, N_USER_RELAZMS, USER_RELAZMS,        & ! Inputs
            FLUX_FACTOR, NBEAMS, BEAM_SZAS, EARTH_RADIUS, HEIGHT_GRID,       & ! Inputs
            DELTAU_INPUT, OMEGA_INPUT, ASYMM_INPUT, D2S_SCALING,             & ! Inputs
            THERMAL_BB_INPUT, LAMBERTIAN_ALBEDO, BRDF_F_0, BRDF_F, UBRDF_F,  & ! Inputs
            EMISSIVITY, SURFBB, SLTERM_ISOTROPIC, SLTERM_F_0,                & ! Inputs  !@@ 2p3 (Add Sleave)
            DO_PROFILE_WFS, DO_SURFACE_WFS, DO_SLEAVE_WFS,                   & ! Inputs  !@@ 2p3 (Add Sleave)
            LAYER_VARY_FLAG, LAYER_VARY_NUMBER, N_SURFACE_WFS, N_SLEAVE_WFS, & ! Inputs
            LSSL_SLTERM_ISOTROPIC, LSSL_SLTERM_F_0,                          & ! Inputs  !@@ 2p3 (New sleave)
            L_DELTAU_INPUT, L_OMEGA_INPUT, L_ASYMM_INPUT, L_D2S_SCALING,     & ! Inputs
            LS_BRDF_F_0, LS_BRDF_F, LS_UBRDF_F, LS_EMISSIVITY,               & ! Inputs
            INTENSITY_TOA, PROFILEWF_TOA, SURFACEWF_TOA,                     & ! In/Out
            INTENSITY_BOA, PROFILEWF_BOA, SURFACEWF_BOA,                     & ! In/Out
            RADLEVEL_UP, RADLEVEL_DN,  N_GEOMETRIES,                         & ! Outputs !@@ 2p2
            PROFJACLEVEL_UP,PROFJACLEVEL_DN,SURFJACLEVEL_UP,SURFJACLEVEL_DN, & ! Outputs !@@ 2p2
            FLUXES_TOA, PROFJACFLUXES_TOA, SURFJACFLUXES_TOA,                & ! Outputs !@@ 2p3
            FLUXES_BOA, PROFJACFLUXES_BOA, SURFJACFLUXES_BOA,                & ! Outputs !@@ 2p3
            STATUS_INPUTCHECK, C_NMESSAGES, C_MESSAGES, C_ACTIONS,           & ! Exception handling
            STATUS_EXECUTION,  E_MESSAGE, E_TRACE_1, E_TRACE_2 )               ! Exception handling

!  Exception handling

          IF ( STATUS_INPUTCHECK .eq. 1 ) THEN
            write(*,'(a,i4)')'INPUT Check failed from Perturbed Profile WF Run # ',T
            write(*,*)' - Number of Messages = ', C_NMESSAGES
            DO k = 1, C_NMESSAGES
              write(*,'(A,I3,A,A)')' - Message # ',K,': ', TRIM(C_MESSAGES(K))
              write(*,'(A,I3,A,A)')' - Action  # ',K,': ', TRIM(C_ACTIONS(K))
            ENDDO
            stop'Test_2S_LPCS program aborted'
          ENDIF
          IF ( STATUS_EXECUTION .eq. 1 ) THEN
            write(*,'(a,i4)')'EXECUTION failed from Perturbed Profile WF Run # ',T
            write(*,*)' - Print 1 Message and 2 Traces'
            write(*,'(A)') TRIM(E_MESSAGE)
            write(*,'(A)') TRIM(E_TRACE_1)
            write(*,'(A)') TRIM(E_TRACE_2)
            stop'Test_2S_LPCS program aborted'
          ENDIF

!  Save the last results produced by the OMP thread

          IF (wvn == (TID+1)*numwvn/OMP_NTHREADS) THEN
            Q = 1
            INTENSITY_TOA_PROF_PT(1:N_GEOMETRIES,I,Q,T,TID+1) = INTENSITY_TOA(1:N_GEOMETRIES)
            INTENSITY_BOA_PROF_PT(1:N_GEOMETRIES,I,Q,T,TID+1) = INTENSITY_BOA(1:N_GEOMETRIES)

            if ( do_2s_LEVELOUT ) THEN
               RADLEVEL_UP_PROF_PT(1:N_GEOMETRIES,0:NLAYERS,I,Q,T,TID+1) = RADLEVEL_UP(1:N_GEOMETRIES,0:NLAYERS)
               RADLEVEL_DN_PROF_PT(1:N_GEOMETRIES,0:NLAYERS,I,Q,T,TID+1) = RADLEVEL_DN(1:N_GEOMETRIES,0:NLAYERS)
            endif

            FLUXES_TOA_PROF_PT(1:NBEAMS,:,I,Q,T,TID+1) = FLUXES_TOA(1:NBEAMS,:)
            FLUXES_BOA_PROF_PT(1:NBEAMS,:,I,Q,T,TID+1) = FLUXES_BOA(1:NBEAMS,:)
          END IF

!  End layer loop

        enddo

!  End thread loop

      enddo


!  FD calculation 3 : Perturbation on aero ext tau
!  ===============================================

      !write(*,*)
      !write(*,*)'Doing OMP thread ',TID,' Doing FD calculation 3: perturb aero ext tau'

!  Start thread loop

      do thread = 1, nthreads
!      do thread = 3, 3
        t = thread

!  Thread 1: 2S, Solar only   , no   delta-M scaling
!  Thread 2: 2S, Solar only   , with delta-M scaling
!  Thread 3: 2S, Thermal only , no   delta-M scaling
!  Thread 4: 2S, Thermal only , with delta-M scaling
!  Thread 5: 2S, Solar+Thermal, no   delta-M scaling
!  Thread 6: 2S, Solar+Thermal, with delta-M scaling

        call thread_setter &
          ( thread, n_user_obsgeoms, & 
            do_solar_sources, do_thermal_emission, &
            do_d2s_scaling, do_surface_emission,   &
            nbeams, n_user_relazms )
        !write(*,*)'Doing thread # ', thread
        !write(*,'(1x,a,i1,a,i1)') 'Doing OMP thread ',TID,' LIDORT thread ',thread

        column = 1.0d0

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
          !emissivity    =  1.0d0 - albedo_save
        else
          surfbb        = 0.0d0
          emissivity    = 0.0d0
        end if

!  Start layer loop

        do i=1,nlayers

          if (i <= n6) then

!  Save the last results produced by the OMP thread

            IF (wvn == (TID+1)*numwvn/OMP_NTHREADS) THEN
              Q = 2

! mick fix 1/16/2015 - set N_GEOMETRIES by hand since TWOSTREAM_LPS_MASTER not called in this case
              N_GEOMETRIES = N_GEOMETRIES_SAVE

              INTENSITY_TOA_PROF_PT(1:N_GEOMETRIES,I,Q,T,TID+1) = INTENSITY_TOA_BAS(1:N_GEOMETRIES,T,TID+1)
              INTENSITY_BOA_PROF_PT(1:N_GEOMETRIES,I,Q,T,TID+1) = INTENSITY_BOA_BAS(1:N_GEOMETRIES,T,TID+1)

              if ( do_2s_LEVELOUT ) THEN
                RADLEVEL_UP_PROF_PT(1:N_GEOMETRIES,0:NLAYERS,I,Q,T,TID+1) = RADLEVEL_UP_BAS(1:N_GEOMETRIES,0:NLAYERS,T,TID+1)
                RADLEVEL_DN_PROF_PT(1:N_GEOMETRIES,0:NLAYERS,I,Q,T,TID+1) = RADLEVEL_DN_BAS(1:N_GEOMETRIES,0:NLAYERS,T,TID+1)
              endif

              FLUXES_TOA_PROF_PT(1:NBEAMS,:,I,Q,T,TID+1) = FLUXES_TOA_BAS(1:NBEAMS,:,T,TID+1)
              FLUXES_BOA_PROF_PT(1:NBEAMS,:,I,Q,T,TID+1) = FLUXES_BOA_BAS(1:NBEAMS,:,T,TID+1)
            END IF

          else

!  Create optical properties, every thread

            do n = 1, n6
              molsca = molomg(n) * molext(n)
              molabs = molext(n) - molsca

              deltau_input(n) = molext(n)
              omega_input(n)  = molomg(n)
              asymm_input(n)  = 0.0d0

              m1     = raymoms(2,n) * molsca ; m2 = 0.0d0
              d2s_scaling(n)  = ( m1 + m2 ) / molsca / 5.0d0
            enddo

            parcel = taer / ( height_grid(n6) - height_grid(nlayers) )
            do n = n6 + 1, nlayers
              molsca = molomg(n) * molext(n)
              molabs = molext(n) - molsca

              !1st try
              if (n == i) then
                aerext = Parcel * ( height_grid(n-1) - height_grid(n) ) * epsfac
              else
                aerext = Parcel * ( height_grid(n-1) - height_grid(n) )
              endif
              aersca = aerext * waer

              !2nd try
              !aerext = Parcel * ( height_grid(n-1) - height_grid(n) )
              !aersca = aerext * waer
              !if (n == i) aerext = Parcel * ( height_grid(n-1) - height_grid(n) ) * epsfac

              totext = molsca + molabs + aerext
              totsca = molsca + aersca
              raywt  = molsca / totsca
              aerwt  = aersca / totsca
              omega  = totsca / totext

              deltau_input(n) = totext
              omega_input(n)  = omega
              asymm_input(n)  = gaer * aerwt

              m1 = raymoms(2,n) * molsca
              m2 = 5.0d0 * gaer * gaer * aersca
              d2s_scaling(n)  = ( m1 + m2 ) / totsca / 5.0d0
            enddo

!  Set up surface BRDF

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

              !Linearized BRDF supplement inputs
              DO_KERNEL_FACTOR_WFS      = .FALSE.
              DO_KERNEL_PARAMS_WFS      = .FALSE.

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
              !      TWOSTREAM_LS_BRDFMASTER using the DO_SHADOW_EFFECT flag

              !Call linearized BRDF supplement
              CALL TWOSTREAM_LS_BRDFMASTER &
            ( MAXBEAMS, MAX_USER_ANGLES, MAX_USER_OBSGEOMS,      & ! Dimensions !@@
              MAXSTREAMS_BRDF, MAX_BRDF_KERNELS,                 & ! Dimensions
              MAX_BRDF_PARAMETERS, MAX_SURFACEWFS,               & ! Dimensions
              DO_SOLAR_SOURCES, DO_USER_OBSGEOMS,                & ! Inputs !@@
              LAMBERTIAN_KERNEL_FLAG,                            & ! Inputs
              DO_SHADOW_EFFECT, DO_SURFACE_EMISSION,             & ! Inputs
              NBEAMS, N_USER_ANGLES, N_USER_OBSGEOMS,            & ! Inputs !@@
              BEAM_SZAS, USER_ANGLES, USER_OBSGEOMS,             & ! Inputs !@@
              STREAM_VALUE, NSTREAMS_BRDF,                       & ! Inputs
              N_BRDF_KERNELS, WHICH_BRDF, BRDF_FACTORS,          & ! Inputs
              N_BRDF_PARAMETERS, BRDF_PARAMETERS,                & ! Inputs
              DO_KERNEL_FACTOR_WFS, DO_KERNEL_PARAMS_WFS,        & ! Inputs
              DO_KPARAMS_DERIVS, N_SURFACE_WFS,                  & ! Outputs
              N_KERNEL_FACTOR_WFS, N_KERNEL_PARAMS_WFS,          & ! Outputs
              BRDF_F_0, BRDF_F, UBRDF_F, EMISSIVITY,             & ! Outputs
              LS_BRDF_F_0, LS_BRDF_F, LS_UBRDF_F, LS_EMISSIVITY, & ! Outputs
              STATUS_BRDFSUP, MESSAGE_BRDF, ACTION_BRDF )          ! Outputs
            else
              !Lambertian surface

              LAMBERTIAN_ALBEDO = albedo_save
              emissivity = 1.0d0 - albedo_save
            end if

!  Exception handling on BRDF

            if ( DO_BRDF_SURFACE ) then
               if ( STATUS_BRDFSUP .ne. 0 ) then
                  write(*,'(a,i4)')'BRDF Lin supplement Check failed from Perturbed Profile WF Run # ',T
                  write(*,*)' - Print 1 Message and 1 Action'
                  write(*,'(A)') TRIM(MESSAGE_BRDF)
                  write(*,'(A)') TRIM(ACTION_BRDF)
                  stop'Test_2S_LPCS program aborted'
               endif
            endif

!  Call to Linearized model

            CALL TWOSTREAM_LPS_MASTER &
            ( MAXLAYERS, MAXTOTAL, MAXMESSAGES, MAXBEAMS, MAX_GEOMETRIES,      & ! Dimensions
              MAX_USER_ANGLES, MAX_USER_RELAZMS, MAX_USER_OBSGEOMS,            & ! Dimensions !@@ 2p1
              MAX_ATMOSWFS, MAX_SURFACEWFS, MAX_SLEAVEWFS,                     & ! Dimensions !@@ 2p3 (Add Sleave)
              DO_UPWELLING, DO_DNWELLING, DO_PLANE_PARALLEL, DO_2S_LEVELOUT,   & ! Inputs     !@@ 2p2
              DO_MVOUT_ONLY, DO_ADDITIONAL_MVOUT,                              & ! Inputs     !@@ 2p3
              DO_SOLAR_SOURCES, DO_THERMAL_EMISSION, DO_SURFACE_EMISSION,      & ! Inputs
              DO_D2S_SCALING, DO_BRDF_SURFACE, DO_USER_OBSGEOMS,               & ! Inputs     !@@ 2p1
              DO_SURFACE_LEAVING, DO_SL_ISOTROPIC, DO_PENTADIAG_INVERSE,       & ! Input !@@ 2p3   6/25/14
              BVPINDEX, BVPSCALEFACTOR, TAYLOR_ORDER, TAYLOR_SMALL,            & ! Input !@@ 2p3/4 6/25/14, 1/7/15
              NLAYERS, NTOTAL, STREAM_VALUE, N_USER_OBSGEOMS, USER_OBSGEOMS,   & ! Inputs     !@@ 2p1
              N_USER_ANGLES, USER_ANGLES, N_USER_RELAZMS, USER_RELAZMS,        & ! Inputs
              FLUX_FACTOR, NBEAMS, BEAM_SZAS, EARTH_RADIUS, HEIGHT_GRID,       & ! Inputs
              DELTAU_INPUT, OMEGA_INPUT, ASYMM_INPUT, D2S_SCALING,             & ! Inputs
              THERMAL_BB_INPUT, LAMBERTIAN_ALBEDO, BRDF_F_0, BRDF_F, UBRDF_F,  & ! Inputs
              EMISSIVITY, SURFBB, SLTERM_ISOTROPIC, SLTERM_F_0,                & ! Inputs  !@@ 2p3 (Add Sleave)
              DO_PROFILE_WFS, DO_SURFACE_WFS, DO_SLEAVE_WFS,                   & ! Inputs  !@@ 2p3 (Add Sleave)
              LAYER_VARY_FLAG, LAYER_VARY_NUMBER, N_SURFACE_WFS, N_SLEAVE_WFS, & ! Inputs
              LSSL_SLTERM_ISOTROPIC, LSSL_SLTERM_F_0,                          & ! Inputs  !@@ 2p3 (New sleave)
              L_DELTAU_INPUT, L_OMEGA_INPUT, L_ASYMM_INPUT, L_D2S_SCALING,     & ! Inputs
              LS_BRDF_F_0, LS_BRDF_F, LS_UBRDF_F, LS_EMISSIVITY,               & ! Inputs
              INTENSITY_TOA, PROFILEWF_TOA, SURFACEWF_TOA,                     & ! In/Out
              INTENSITY_BOA, PROFILEWF_BOA, SURFACEWF_BOA,                     & ! In/Out
              RADLEVEL_UP, RADLEVEL_DN,  N_GEOMETRIES,                         & ! Outputs !@@ 2p2
              PROFJACLEVEL_UP,PROFJACLEVEL_DN,SURFJACLEVEL_UP,SURFJACLEVEL_DN, & ! Outputs !@@ 2p2
              FLUXES_TOA, PROFJACFLUXES_TOA, SURFJACFLUXES_TOA,                & ! Outputs !@@ 2p3
              FLUXES_BOA, PROFJACFLUXES_BOA, SURFJACFLUXES_BOA,                & ! Outputs !@@ 2p3
              STATUS_INPUTCHECK, C_NMESSAGES, C_MESSAGES, C_ACTIONS,           & ! Exception handling
              STATUS_EXECUTION,  E_MESSAGE, E_TRACE_1, E_TRACE_2 )               ! Exception handling

!  Exception handling

            IF ( STATUS_INPUTCHECK .eq. 1 ) THEN
              write(*,'(a,i4)')'INPUT Check failed from Perturbed Profile WF Run # ',T
              write(*,*)' - Number of Messages = ', C_NMESSAGES
              DO k = 1, C_NMESSAGES
                write(*,'(A,I3,A,A)')' - Message # ',K,': ', TRIM(C_MESSAGES(K))
                write(*,'(A,I3,A,A)')' - Action  # ',K,': ', TRIM(C_ACTIONS(K))
              ENDDO
              stop'Test_2S_LPCS program aborted'
            ENDIF
            IF ( STATUS_EXECUTION .eq. 1 ) THEN
              write(*,'(a,i4)')'EXECUTION failed from Perturbed Profile WF Run # ',T
              write(*,*)' - Print 1 Message and 2 Traces'
              write(*,'(A)') TRIM(E_MESSAGE)
              write(*,'(A)') TRIM(E_TRACE_1)
              write(*,'(A)') TRIM(E_TRACE_2)
              stop'Test_2S_LPCS program aborted'
            ENDIF

!  Save the last results produced by the OMP thread

            IF (wvn == (TID+1)*numwvn/OMP_NTHREADS) THEN
              Q = 2
              INTENSITY_TOA_PROF_PT(1:N_GEOMETRIES,I,Q,T,TID+1) = INTENSITY_TOA(1:N_GEOMETRIES)
              INTENSITY_BOA_PROF_PT(1:N_GEOMETRIES,I,Q,T,TID+1) = INTENSITY_BOA(1:N_GEOMETRIES)

              if ( do_2s_LEVELOUT ) THEN
                 RADLEVEL_UP_PROF_PT(1:N_GEOMETRIES,0:NLAYERS,I,Q,T,TID+1) = RADLEVEL_UP(1:N_GEOMETRIES,0:NLAYERS)
                 RADLEVEL_DN_PROF_PT(1:N_GEOMETRIES,0:NLAYERS,I,Q,T,TID+1) = RADLEVEL_DN(1:N_GEOMETRIES,0:NLAYERS)
              endif

              FLUXES_TOA_PROF_PT(1:NBEAMS,:,I,Q,T,TID+1) = FLUXES_TOA(1:NBEAMS,:)
              FLUXES_BOA_PROF_PT(1:NBEAMS,:,I,Q,T,TID+1) = FLUXES_BOA(1:NBEAMS,:)
            END IF

          endif

!  End layer loop

        enddo

!  End thread loop

      enddo

!  End "wavenumber" loop

      enddo

!$OMP END DO

      !IF (TID == 0) call cpu_time(omp_e2)

!  End parallel region

!$OMP END PARALLEL

!  Timing tests (OpenMP)
!  =====================

      call cpu_time(omp_e2)

      if (omp_nthreads <= n_core) then
        time_divider = omp_nthreads
      else
        time_divider = n_core
      endif

      write(*,*)
      write(*,'(15x,a)')                    'Timing report'
      write(*,'(4(1x,a))') 'Numwvn','# OMP Threads','Time (sec)'
      write(*,'(4(1x,a))') '------','-------------','----------'
      write(*,'(1x,i5,8x,i1,7x,f6.2)') &
        numwvn, omp_nthreads, (omp_e2 - omp_e1)/real(time_divider)

!  Write Part 1 Results
!  --------------------

      CALL WRITE_2S_LPS_OUTPUT()

!     stop'profile only'

656   continue

      write(*,*)
      write(*,*)'**************************************'
      write(*,*)

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!                   Part 2: COLUMN/SURFACE WFS TESTS                  @
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  Initialize saved output for Intensity and associated Jacobians

      INTENSITY_TOA_BAS     = 0.0d0
      INTENSITY_TOA_COL_PT  = 0.0d0
      INTENSITY_TOA_SURF_PT = 0.0d0
      COLUMNWF_TOA_BAS      = 0.0d0
      SURFACEWF_TOA_BAS     = 0.0d0

      INTENSITY_BOA_BAS     = 0.0d0
      INTENSITY_BOA_COL_PT  = 0.0d0
      INTENSITY_BOA_SURF_PT = 0.0d0
      COLUMNWF_BOA_BAS      = 0.0d0
      SURFACEWF_BOA_BAS     = 0.0d0

!  Initialize saved output for Level Radiances and associated Jacobians

      RADLEVEL_UP_BAS       = 0.0d0
      RADLEVEL_UP_COL_PT    = 0.0d0
      RADLEVEL_UP_SURF_PT   = 0.0d0
      COLJACLEVEL_UP_BAS    = 0.0d0
      SURFJACLEVEL_UP_BAS   = 0.0d0

      RADLEVEL_DN_BAS       = 0.0d0
      RADLEVEL_DN_COL_PT    = 0.0d0
      RADLEVEL_DN_SURF_PT   = 0.0d0
      COLJACLEVEL_DN_BAS    = 0.0d0
      SURFJACLEVEL_DN_BAS   = 0.0d0

!  Initialize saved output for Fluxes and associated Jacobians

      FLUXES_TOA_BAS        = 0.0d0
      FLUXES_TOA_COL_PT     = 0.0d0
      FLUXES_TOA_SURF_PT    = 0.0d0
      COLJACFLUXES_TOA_BAS  = 0.0d0
      SURFJACFLUXES_TOA_BAS = 0.0d0

      FLUXES_BOA_BAS        = 0.0d0
      FLUXES_BOA_COL_PT     = 0.0d0
      FLUXES_BOA_SURF_PT    = 0.0d0
      COLJACFLUXES_BOA_BAS  = 0.0d0
      SURFJACFLUXES_BOA_BAS = 0.0d0

!     Begin parallel region

!$OMP PARALLEL IF(OMP_MAXTHREADS > 1) &
!$OMP   DEFAULT (PRIVATE) &
!$OMP   SHARED (OMP_NTHREADS, &
!$OMP     DO_UPWELLING, DO_DNWELLING, &
!$OMP     DO_SURFACE_LEAVING, DO_SL_ISOTROPIC, DO_SLEAVE_WFS, N_SLEAVE_WFS, &
!$OMP     BVPINDEX, BVPSCALEFACTOR, TAYLOR_ORDER, TAYLOR_SMALL, &
!$OMP     NUMWVN, NTHREADS, NLAYERS, NTOTAL, &
!$OMP     N_USER_OBSGEOMS, N_USER_ANGLES, NSTREAMS_BRDF, &
!$OMP     BEAM_SZAS, USER_ANGLES, USER_RELAZMS, USER_OBSGEOMS, &
!$OMP     STREAM_VALUE, FLUX_FACTOR, EARTH_RADIUS, &
!$OMP     HEIGHT_GRID, MOLEXT, MOLOMG, RAYMOMS, N6, GAER, WAER, TAER, &
!$OMP     NTEMPS, AIRTEMPS, SURFTEMP, &
!$OMP     WINDSPEED_SAVE, ALBEDO_SAVE, BRDF_PAR, eps, &
!$OMP     INTENSITY_TOA_BAS, INTENSITY_BOA_BAS, &
!$OMP     COLUMNWF_TOA_BAS,  COLUMNWF_BOA_BAS,  INTENSITY_TOA_COL_PT,  INTENSITY_BOA_COL_PT, &
!$OMP     SURFACEWF_TOA_BAS, SURFACEWF_BOA_BAS, INTENSITY_TOA_SURF_PT, INTENSITY_BOA_SURF_PT, &
!$OMP     RADLEVEL_UP_BAS, RADLEVEL_DN_BAS, & 
!$OMP     COLJACLEVEL_UP_BAS,  COLJACLEVEL_DN_BAS,  RADLEVEL_UP_COL_PT,  RADLEVEL_DN_COL_PT, &
!$OMP     SURFJACLEVEL_UP_BAS, SURFJACLEVEL_DN_BAS, RADLEVEL_UP_SURF_PT, RADLEVEL_DN_SURF_PT, &
!$OMP     FLUXES_TOA_BAS, FLUXES_BOA_BAS, &
!$OMP     COLJACFLUXES_TOA_BAS,  COLJACFLUXES_BOA_BAS,  FLUXES_TOA_COL_PT,  FLUXES_BOA_COL_PT, &
!$OMP     SURFJACFLUXES_TOA_BAS, SURFJACFLUXES_BOA_BAS, FLUXES_TOA_SURF_PT, FLUXES_BOA_SURF_PT, &
!$OMP     N_GEOMETRIES_SAVE) &
!$OMP   FIRSTPRIVATE (DO_SOLAR_SOURCES, DO_THERMAL_EMISSION, DO_SURFACE_EMISSION, &
!$OMP     DO_BRDF_SURFACE, DO_FULLQUADRATURE, DO_D2S_SCALING, &
!$OMP     DO_PLANE_PARALLEL, DO_USER_OBSGEOMS, DO_2S_LEVELOUT, &
!$OMP     DO_MVOUT_ONLY, DO_ADDITIONAL_MVOUT, DO_PENTADIAG_INVERSE, &
!$OMP     NBEAMS, N_USER_RELAZMS)

!     Obtain thread number

      TID = OMP_GET_THREAD_NUM()
      !TID = 0
      
!     Obtain and display total number of threads and local thread

      IF (TID == 0) THEN
        OMP_NTHREADS = OMP_GET_NUM_THREADS()
        !OMP_NTHREADS = 1
        !write(*,*)
        !write(*,'(1x,a,i1)') 'Total number of threads (inside parallel region) = ', OMP_NTHREADS
        write(*,'(1x,a,i1,a)') 'Running 2S driver with ',OMP_NTHREADS,' thread(s)'
        write(*,*)
      END IF

      write(*,'(1x,a,i1,a)') 'Thread ',TID,' beginning its wavenumber loop'

!$OMP DO

!  Begin "wavenumber" loop

      do wvn = 1, numwvn

      !write(*,'(1x,a,i1,a,i5)') 'Thread  = ',TID,' wvn = ',wvn
      write(*,*)
      write(*,*)'Doing OMP thread ',TID,' Column WF baseline: Intensities, 2 Column WFs + 1 Surface WF'
      write(*,*) 'wvn = ',wvn

!  Initialize some inputs

      deltau_input = 0.0d0
      omega_input  = 0.0d0
      asymm_input  = 0.0d0
      d2s_scaling  = 0.0d0

      lambertian_albedo = 0.0d0

      BRDF_F_0 = 0.0d0
      BRDF_F   = 0.0d0
      UBRDF_F  = 0.0d0

      L_deltau_input = 0.0d0
      L_omega_input  = 0.0d0
      L_asymm_input  = 0.0d0
      L_d2s_scaling  = 0.0d0

      LS_BRDF_F_0 = 0.0d0
      LS_BRDF_F   = 0.0d0
      LS_UBRDF_F  = 0.0d0

!  Control for 2 column WFs

      do_profile_wfs = .false.
      do_column_wfs  = .true.
      do_surface_wfs = .true.
      do n = 1, nlayers
         layer_vary_flag(n)   = .false.
         layer_vary_number(n) = 0
      enddo
      n_column_wfs  = 2
      n_surface_wfs = 1

!  Control for SLEAVE

      DO_SURFACE_LEAVING = .FALSE.
      DO_SL_ISOTROPIC    = .FALSE.
      SLTERM_ISOTROPIC   = 0.0D0
      SLTERM_F_0         = 0.0D0
      DO_SLEAVE_WFS      = .FALSE.
      N_SLEAVE_WFS       = 0
      LSSL_SLTERM_ISOTROPIC = 0.0D0
      LSSL_SLTERM_F_0       = 0.0D0

!  Baseline calculation 2 : RADIANCE + COLUMN/SURFACE WFS
!  ======================================================

!  Start thread loop

      do thread = 1, nthreads
!      do thread = 1, 1
        t = thread

!  Thread 1: 2S, Solar only   , no   delta-M scaling
!  Thread 2: 2S, Solar only   , with delta-M scaling
!  Thread 3: 2S, Thermal only , no   delta-M scaling
!  Thread 4: 2S, Thermal only , with delta-M scaling
!  Thread 5: 2S, Solar+Thermal, no   delta-M scaling
!  Thread 6: 2S, Solar+Thermal, with delta-M scaling

        call thread_setter &
          ( thread, n_user_obsgeoms, & 
            do_solar_sources, do_thermal_emission, &
            do_d2s_scaling, do_surface_emission,   &
            nbeams, n_user_relazms )
        !write(*,*)'Doing thread # ', thread
        !write(*,'(1x,a,i1,a,i1)') 'Doing OMP thread ',TID,' LIDORT thread ',thread

        column = 1.0d0

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
          !emissivity    =  1.0d0 - albedo_save
          !ls_emissivity = -1.0d0
        else
          surfbb        = 0.0d0
          emissivity    = 0.0d0
          ls_emissivity = 0.0d0
        end if

!  Create optical properties, every thread

        do n = 1, n6
          deltau_input(n) = molext(n)
          omega_input(n)  = molomg(n)
          asymm_input(n)  = 0.0d0

          molsca = molomg(n) * molext(n)
          molabs = molext(n) - molsca
          m1     = raymoms(2,n) * molsca ; m2 = 0.0d0
          d2s_scaling(n)  = ( m1 + m2 ) / molsca / 5.0d0


          !1st linearized oprop set (partially normalized) -
          !  par=1 wrt mol abs tau
          omega = omega_input(n)
          l_deltau_input(n,1) = molabs
          l_omega_input (n,1) = - molabs * omega/deltau_input(n)
        enddo

        parcel = taer / ( height_grid(n6) - height_grid(nlayers) )
        do n = n6 + 1, nlayers
          aerext = Parcel * ( height_grid(n-1) - height_grid(n) )
          aersca = aerext * waer
          molsca = molomg(n) * molext(n)
          molabs = molext(n) - molsca

          totext = molext(n) + aerext
          totsca = molsca    + aersca
          raywt  = molsca / totsca
          aerwt  = aersca / totsca
          omega  = totsca / totext

          deltau_input(n) = totext
          omega_input(n)  = omega
          asymm_input(n)  = gaer * aersca / totsca

          m1 = raymoms(2,n) * molsca
          m2 = 5.0d0 * gaer * gaer * aersca
          d2s_scaling(n)  = ( m1 + m2 ) / totsca / 5.0d0


          !1st linearized oprop set (partially normalized) -
          !  par=1 wrt mol abs tau
          l_deltau_input(n,1) = molabs
          l_omega_input (n,1) = - molabs * omega /deltau_input(n)


          !2nd linearized oprop set (partially normalized) -
          !  par=2 wrt aero ext tau


! @@@@@@@@@@@@@@@@@@@@@@
!  This set was wrong
!          l_deltau_input(n,2) = aerext
!          l_omega_input(n,2)  = aerext * (1.0d0 - omega) / deltau_input(n)
!          l_asymm_input(n,2)  = aerext * (gaer - asymm_input(n)) / totsca
!          l_d2s_scaling(n,2)  = aerext * (gaer*gaer - d2s_scaling(n)) / totsca

!  Rob fixed, 14 nov 12
          l_deltau_input(n,2) = aerext
          l_omega_input(n,2)  = aerext * (waer - omega) / deltau_input(n)
          l_asymm_input(n,2)  = aerwt  * (gaer - asymm_input(n))
          l_d2s_scaling(n,2)  = aerwt  * (gaer*gaer - d2s_scaling(n))
! @@@@@@@@@@@@@@@@@@@@@@

        enddo

!  Set up surface BRDF

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

          !Linearized BRDF supplement inputs
          DO_KERNEL_FACTOR_WFS      = .FALSE.
          DO_KERNEL_PARAMS_WFS      = .FALSE.

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
          !      TWOSTREAM_LS_BRDFMASTER using the DO_SHADOW_EFFECT flag

          !Linearized BRDF supplement inputs
          DO_KERNEL_PARAMS_WFS(1,1) = .TRUE.

          !Save BRDF slope parameter for perturbation comparison below
          BRDF_PAR = BRDF_PARAMETERS(1,1)

          !Call linearized BRDF supplement
          CALL TWOSTREAM_LS_BRDFMASTER &
            ( MAXBEAMS, MAX_USER_ANGLES, MAX_USER_OBSGEOMS,      & ! Dimensions !@@
              MAXSTREAMS_BRDF, MAX_BRDF_KERNELS,                 & ! Dimensions
              MAX_BRDF_PARAMETERS, MAX_SURFACEWFS,               & ! Dimensions
              DO_SOLAR_SOURCES, DO_USER_OBSGEOMS,                & ! Inputs !@@
              LAMBERTIAN_KERNEL_FLAG,                            & ! Inputs
              DO_SHADOW_EFFECT, DO_SURFACE_EMISSION,             & ! Inputs
              NBEAMS, N_USER_ANGLES, N_USER_OBSGEOMS,            & ! Inputs !@@
              BEAM_SZAS, USER_ANGLES, USER_OBSGEOMS,             & ! Inputs !@@
              STREAM_VALUE, NSTREAMS_BRDF,                       & ! Inputs
              N_BRDF_KERNELS, WHICH_BRDF, BRDF_FACTORS,          & ! Inputs
              N_BRDF_PARAMETERS, BRDF_PARAMETERS,                & ! Inputs
              DO_KERNEL_FACTOR_WFS, DO_KERNEL_PARAMS_WFS,        & ! Inputs
              DO_KPARAMS_DERIVS, N_SURFACE_WFS,                  & ! Outputs
              N_KERNEL_FACTOR_WFS, N_KERNEL_PARAMS_WFS,          & ! Outputs
              BRDF_F_0, BRDF_F, UBRDF_F, EMISSIVITY,             & ! Outputs
              LS_BRDF_F_0, LS_BRDF_F, LS_UBRDF_F, LS_EMISSIVITY, & ! Outputs
              STATUS_BRDFSUP, MESSAGE_BRDF, ACTION_BRDF )          ! Outputs

        else
          !Lambertian surface

          BRDF_PAR = albedo_save
          LAMBERTIAN_ALBEDO = albedo_save
          emissivity    =  1.0d0 - albedo_save
          ls_emissivity = -1.0d0
        end if

!  Exception handling on BRDF

        if ( DO_BRDF_SURFACE ) then
           if ( STATUS_BRDFSUP .ne. 0 ) then
              write(*,'(a,i4)')'BRDF Lin supplement Check failed from Baseline Column WF Run # ',T
              write(*,*)' - Print 1 Message and 1 Action'
              write(*,'(A)') TRIM(MESSAGE_BRDF)
              write(*,'(A)') TRIM(ACTION_BRDF)
              stop'Test_2S_LPCS program aborted'
           endif
        endif

!  Call to Linearized model

        CALL TWOSTREAM_LCS_MASTER &
        ( MAXLAYERS, MAXTOTAL, MAXMESSAGES, MAXBEAMS, MAX_GEOMETRIES,      & ! Dimensions
          MAX_USER_ANGLES, MAX_USER_RELAZMS, MAX_USER_OBSGEOMS,            & ! Dimensions !@@ 2p1
          MAX_ATMOSWFS, MAX_SURFACEWFS, MAX_SLEAVEWFS,                     & ! Dimensions !@@ 2p3 (Add Sleave)
          DO_UPWELLING, DO_DNWELLING, DO_PLANE_PARALLEL, DO_2S_LEVELOUT,   & ! Inputs     !@@ 2p2
          DO_MVOUT_ONLY, DO_ADDITIONAL_MVOUT,                              & ! Inputs     !@@ 2p3
          DO_SOLAR_SOURCES, DO_THERMAL_EMISSION, DO_SURFACE_EMISSION,      & ! Inputs
          DO_D2S_SCALING, DO_BRDF_SURFACE, DO_USER_OBSGEOMS,               & ! Inputs     !@@ 2p1
          DO_SURFACE_LEAVING, DO_SL_ISOTROPIC, DO_PENTADIAG_INVERSE,       & ! Input !@@ 2p3   6/25/14
          BVPINDEX, BVPSCALEFACTOR, TAYLOR_ORDER, TAYLOR_SMALL,            & ! Input !@@ 2p3/4 6/25/14, 1/7/15
          NLAYERS, NTOTAL, STREAM_VALUE, N_USER_OBSGEOMS, USER_OBSGEOMS,   & ! Inputs     !@@ 2p1
          N_USER_ANGLES, USER_ANGLES, N_USER_RELAZMS, USER_RELAZMS,        & ! Inputs
          FLUX_FACTOR, NBEAMS, BEAM_SZAS, EARTH_RADIUS, HEIGHT_GRID,       & ! Inputs
          DELTAU_INPUT, OMEGA_INPUT, ASYMM_INPUT, D2S_SCALING,             & ! Inputs
          THERMAL_BB_INPUT, LAMBERTIAN_ALBEDO, BRDF_F_0, BRDF_F, UBRDF_F,  & ! Inputs
          EMISSIVITY, SURFBB, SLTERM_ISOTROPIC, SLTERM_F_0,                & ! Inputs  !@@ 2p3 (Add Sleave)
          DO_COLUMN_WFS, DO_SURFACE_WFS, DO_SLEAVE_WFS,                    & ! Inputs  !@@ 2p3 (Add Sleave)
          N_COLUMN_WFS, N_SURFACE_WFS, N_SLEAVE_WFS,                       & ! Inputs  !@@ 2p3 (Add Sleave)
          LSSL_SLTERM_ISOTROPIC, LSSL_SLTERM_F_0,                          & ! Inputs  !@@ 2p3 (New sleave)
          L_DELTAU_INPUT, L_OMEGA_INPUT, L_ASYMM_INPUT, L_D2S_SCALING,     & ! Inputs
          LS_BRDF_F_0, LS_BRDF_F, LS_UBRDF_F, LS_EMISSIVITY,               & ! Inputs
          INTENSITY_TOA, COLUMNWF_TOA, SURFACEWF_TOA,                      & ! Outputs
          INTENSITY_BOA, COLUMNWF_BOA, SURFACEWF_BOA,                      & ! Outputs
          RADLEVEL_UP, RADLEVEL_DN,  N_GEOMETRIES,                         & ! Outputs !@@ 2p2
          COLJACLEVEL_UP, COLJACLEVEL_DN, SURFJACLEVEL_UP, SURFJACLEVEL_DN,& ! Outputs !@@ 2p2
          FLUXES_TOA, COLJACFLUXES_TOA, SURFJACFLUXES_TOA,                 & ! Outputs !@@ 2p3
          FLUXES_BOA, COLJACFLUXES_BOA, SURFJACFLUXES_BOA,                 & ! Outputs !@@ 2p3
          STATUS_INPUTCHECK, C_NMESSAGES, C_MESSAGES, C_ACTIONS,           & ! Exception handling
          STATUS_EXECUTION,  E_MESSAGE, E_TRACE_1, E_TRACE_2 )               ! Exception handling

!  Exception handling

        IF ( STATUS_INPUTCHECK .eq. 1 ) THEN
          write(*,'(a,i4)')'INPUT Check failed from Baseline Column WF Run # ',T
          write(*,*)' - Number of Messages = ', C_NMESSAGES
          DO k = 1, C_NMESSAGES
            write(*,'(A,I3,A,A)')' - Message # ',K,': ', TRIM(C_MESSAGES(K))
            write(*,'(A,I3,A,A)')' - Action  # ',K,': ', TRIM(C_ACTIONS(K))
          ENDDO
          stop'Test_2S_LPCS program aborted'
        ENDIF
        IF ( STATUS_EXECUTION .eq. 1 ) THEN
          write(*,'(a,i4)')'EXECUTION failed from Baseline Column WF Run # ',T
          write(*,*)' - Print 1 Message and 2 Traces'
          write(*,'(A)') TRIM(E_MESSAGE)
          write(*,'(A)') TRIM(E_TRACE_1)
          write(*,'(A)') TRIM(E_TRACE_2)
          stop'Test_2S_LPCS program aborted'
        ENDIF

!  Save the last results produced by the OMP thread

        IF (wvn == (TID+1)*numwvn/OMP_NTHREADS) THEN
          INTENSITY_TOA_BAS ( 1:N_GEOMETRIES, T, TID+1 ) = &
              INTENSITY_TOA ( 1:N_GEOMETRIES )
          INTENSITY_BOA_BAS ( 1:N_GEOMETRIES, T, TID+1 ) = &
              INTENSITY_BOA ( 1:N_GEOMETRIES )

          COLUMNWF_TOA_BAS ( 1:N_GEOMETRIES, :, T, TID+1 ) = &
              COLUMNWF_TOA ( 1:N_GEOMETRIES, : )
          COLUMNWF_BOA_BAS ( 1:N_GEOMETRIES, :, T, TID+1 ) = &
              COLUMNWF_BOA ( 1:N_GEOMETRIES, : )

          SURFACEWF_TOA_BAS ( 1:N_GEOMETRIES, 1:N_SURFACE_WFS, T, TID+1 ) = &
              SURFACEWF_TOA ( 1:N_GEOMETRIES, 1:N_SURFACE_WFS )
          SURFACEWF_BOA_BAS ( 1:N_GEOMETRIES, 1:N_SURFACE_WFS, T, TID+1 ) = &
              SURFACEWF_BOA ( 1:N_GEOMETRIES, 1:N_SURFACE_WFS )

          IF ( do_2s_LEVELOUT ) THEN
            RADLEVEL_UP_BAS ( 1:N_GEOMETRIES, 0:NLAYERS, T, TID+1 ) = &
                RADLEVEL_UP ( 1:N_GEOMETRIES, 0:NLAYERS )
            RADLEVEL_DN_BAS ( 1:N_GEOMETRIES, 0:NLAYERS, T, TID+1 ) = &
                RADLEVEL_DN ( 1:N_GEOMETRIES, 0:NLAYERS )

            COLJACLEVEL_UP_BAS ( 1:N_GEOMETRIES, 0:NLAYERS, :, T, TID+1 ) = &
                COLJACLEVEL_UP ( 1:N_GEOMETRIES, 0:NLAYERS, : )
            COLJACLEVEL_DN_BAS ( 1:N_GEOMETRIES, 0:NLAYERS, :, T, TID+1 ) = &
                COLJACLEVEL_DN ( 1:N_GEOMETRIES, 0:NLAYERS, : )

            SURFJACLEVEL_UP_BAS ( 1:N_GEOMETRIES, 0:NLAYERS, 1:N_SURFACE_WFS, T, TID+1 ) = &
                SURFJACLEVEL_UP ( 1:N_GEOMETRIES, 0:NLAYERS, 1:N_SURFACE_WFS )
            SURFJACLEVEL_DN_BAS ( 1:N_GEOMETRIES, 0:NLAYERS, 1:N_SURFACE_WFS, T, TID+1 ) = &
                SURFJACLEVEL_DN ( 1:N_GEOMETRIES, 0:NLAYERS, 1:N_SURFACE_WFS )
          ENDIF

          FLUXES_TOA_BAS ( 1:NBEAMS, :, T, TID+1 ) = &
              FLUXES_TOA ( 1:NBEAMS, : )
          FLUXES_BOA_BAS ( 1:NBEAMS, :, T, TID+1 ) = &
              FLUXES_BOA ( 1:NBEAMS, : )

          COLJACFLUXES_TOA_BAS ( 1:NBEAMS, :, :, T, TID+1 ) = &
              COLJACFLUXES_TOA ( 1:NBEAMS, :, : )
          COLJACFLUXES_BOA_BAS ( 1:NBEAMS, :, :, T, TID+1 ) = &
              COLJACFLUXES_BOA ( 1:NBEAMS, :, : )

          SURFJACFLUXES_TOA_BAS ( 1:NBEAMS, :, 1:N_SURFACE_WFS, T, TID+1 ) = &
              SURFJACFLUXES_TOA ( 1:NBEAMS, :, 1:N_SURFACE_WFS )
          SURFJACFLUXES_BOA_BAS ( 1:NBEAMS, :, 1:N_SURFACE_WFS, T, TID+1 ) = &
              SURFJACFLUXES_BOA ( 1:NBEAMS, :, 1:N_SURFACE_WFS )
        END IF

!  End thread loop

      enddo


!  Preparing for FD Column Jacobian calculations
!  =============================================

!  Zeroing

      L_deltau_input = 0.0d0
      L_omega_input  = 0.0d0
      L_asymm_input  = 0.0d0
      L_d2s_scaling  = 0.0d0

      LS_BRDF_F_0    = 0.0d0
      LS_BRDF_F      = 0.0d0
      LS_UBRDF_F     = 0.0d0
      LS_EMISSIVITY  = 0.0d0

!  Control for no WFs

      do_profile_wfs = .false.
      do_column_wfs  = .false.
      do_surface_wfs = .false.
      n_column_wfs   = 0
      n_surface_wfs  = 0
      do n = 1, nlayers
         layer_vary_flag(n)   = .false.
         layer_vary_number(n) = 0
      enddo

!  FD perturbation

      !eps = 1.0d-03
      epsfac = 1.0d0 + eps


!  FD calculation 1 : Perturbation on Cox-Munk slope parameter
!  ===========================================================

      !write(*,*)
      !write(*,*)'Doing OMP thread ',TID,' Doing FD calculation 1: perturb Cox-Munk slope parameter'

!  Start thread loop

      do thread = 1, nthreads
!      do thread = 3, 3
        t = thread

!  Thread 1: 2S, Solar only   , no   delta-M scaling
!  Thread 2: 2S, Solar only   , with delta-M scaling
!  Thread 3: 2S, Thermal only , no   delta-M scaling
!  Thread 4: 2S, Thermal only , with delta-M scaling
!  Thread 5: 2S, Solar+Thermal, no   delta-M scaling
!  Thread 6: 2S, Solar+Thermal, with delta-M scaling

        call thread_setter &
          ( thread, n_user_obsgeoms, & 
            do_solar_sources, do_thermal_emission, &
            do_d2s_scaling, do_surface_emission,   &
            nbeams, n_user_relazms )
        !write(*,*)'Doing thread # ', thread
        !write(*,'(1x,a,i1,a,i1)') 'Doing OMP thread ',TID,' LIDORT thread ',thread

        column = 1.0d0

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
          !emissivity    =  1.0d0 - albedo_save
        else
          surfbb        = 0.0d0
          emissivity    = 0.0d0
        end if

!  Create optical properties, every thread

        do n = 1, n6
          deltau_input(n) = molext(n)
          omega_input(n)  = molomg(n)
          asymm_input(n)  = 0.0d0

          molsca = molomg(n) * molext(n)
          molabs = molext(n) - molsca
          m1     = raymoms(2,n) * molsca ; m2 = 0.0d0
          d2s_scaling(n)  = ( m1 + m2 ) / molsca / 5.0d0
        enddo

        parcel = taer / ( height_grid(n6) - height_grid(nlayers) )
        do n = n6 + 1, nlayers
          aerext = Parcel * ( height_grid(n-1) - height_grid(n) )
          molsca = molomg(n) * molext(n)
          molabs = molext(n) - molsca

!          molabs = molabs * 1.001d0
!          aerext = aerext * 1.001d0

          aersca = aerext * waer
          totext = molsca + molabs + aerext
          totsca = molsca + aersca
          raywt  = molsca / totsca
          aerwt  = aersca / totsca
          omega  = totsca / totext

          deltau_input(n) = totext
          omega_input(n)  = omega
          asymm_input(n)  = gaer * aerwt

          m1 = raymoms(2,n) * molsca
          m2 = 5.0d0 * gaer * gaer * aersca
          d2s_scaling(n)  = ( m1 + m2 ) / totsca / 5.0d0
        enddo

!  Set up surface BRDF

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

          !Linearized BRDF supplement inputs
          DO_KERNEL_FACTOR_WFS      = .FALSE.
          DO_KERNEL_PARAMS_WFS      = .FALSE.

          !(2) Now, for the Cox-Munk BRDF, define:

          !Basic BRDF supplement inputs
          N_BRDF_KERNELS            = 1

          WHICH_BRDF(1)             = 9 !Cox-Munk
          BRDF_FACTORS(1)           = 1.0D0
          N_BRDF_PARAMETERS(1)      = 3
          WINDSPEED(T)              = windspeed_save
          BRDF_PARAMETERS(1,1)      = (0.003d0+0.00512d0*WINDSPEED(T))*epsfac !Slope squared
          BRDF_PARAMETERS(1,2)      = 1.33D0*1.33D0 !Square of refractive index
          !Note: The Cox-Munk shadow effect parameter (the 3rd
          !      parameter BRDF_PARAMETERS(1,3)) is set internally by
          !      TWOSTREAM_LS_BRDFMASTER using the DO_SHADOW_EFFECT flag

          !Call linearized BRDF supplement
          CALL TWOSTREAM_LS_BRDFMASTER &
            ( MAXBEAMS, MAX_USER_ANGLES, MAX_USER_OBSGEOMS,      & ! Dimensions !@@
              MAXSTREAMS_BRDF, MAX_BRDF_KERNELS,                 & ! Dimensions
              MAX_BRDF_PARAMETERS, MAX_SURFACEWFS,               & ! Dimensions
              DO_SOLAR_SOURCES, DO_USER_OBSGEOMS,                & ! Inputs !@@
              LAMBERTIAN_KERNEL_FLAG,                            & ! Inputs
              DO_SHADOW_EFFECT, DO_SURFACE_EMISSION,             & ! Inputs
              NBEAMS, N_USER_ANGLES, N_USER_OBSGEOMS,            & ! Inputs !@@
              BEAM_SZAS, USER_ANGLES, USER_OBSGEOMS,             & ! Inputs !@@
              STREAM_VALUE, NSTREAMS_BRDF,                       & ! Inputs
              N_BRDF_KERNELS, WHICH_BRDF, BRDF_FACTORS,          & ! Inputs
              N_BRDF_PARAMETERS, BRDF_PARAMETERS,                & ! Inputs
              DO_KERNEL_FACTOR_WFS, DO_KERNEL_PARAMS_WFS,        & ! Inputs
              DO_KPARAMS_DERIVS, N_SURFACE_WFS,                  & ! Outputs
              N_KERNEL_FACTOR_WFS, N_KERNEL_PARAMS_WFS,          & ! Outputs
              BRDF_F_0, BRDF_F, UBRDF_F, EMISSIVITY,             & ! Outputs
              LS_BRDF_F_0, LS_BRDF_F, LS_UBRDF_F, LS_EMISSIVITY, & ! Outputs
              STATUS_BRDFSUP, MESSAGE_BRDF, ACTION_BRDF )          ! Outputs

        else
          !Lambertian surface

          LAMBERTIAN_ALBEDO = albedo_save * epsfac
          emissivity =  1.0d0 - LAMBERTIAN_ALBEDO
        end if

!  Exception handling on BRDF

        if ( DO_BRDF_SURFACE ) then
           if ( STATUS_BRDFSUP .ne. 0 ) then
              write(*,'(a,i4)')'BRDF Lin supplement Check failed from Perturbed Column WF Run # ',T
              write(*,*)' - Print 1 Message and 1 Action'
              write(*,'(A)') TRIM(MESSAGE_BRDF)
              write(*,'(A)') TRIM(ACTION_BRDF)
              stop'Test_2S_LPCS program aborted'
           endif
        endif

!  Call to Linearized model

        CALL TWOSTREAM_LCS_MASTER &
        ( MAXLAYERS, MAXTOTAL, MAXMESSAGES, MAXBEAMS, MAX_GEOMETRIES,      & ! Dimensions
          MAX_USER_ANGLES, MAX_USER_RELAZMS, MAX_USER_OBSGEOMS,            & ! Dimensions !@@ 2p1
          MAX_ATMOSWFS, MAX_SURFACEWFS, MAX_SLEAVEWFS,                     & ! Dimensions !@@ 2p3 (Add Sleave)
          DO_UPWELLING, DO_DNWELLING, DO_PLANE_PARALLEL, DO_2S_LEVELOUT,   & ! Inputs     !@@ 2p2
          DO_MVOUT_ONLY, DO_ADDITIONAL_MVOUT,                              & ! Inputs     !@@ 2p3
          DO_SOLAR_SOURCES, DO_THERMAL_EMISSION, DO_SURFACE_EMISSION,      & ! Inputs
          DO_D2S_SCALING, DO_BRDF_SURFACE, DO_USER_OBSGEOMS,               & ! Inputs     !@@ 2p1
          DO_SURFACE_LEAVING, DO_SL_ISOTROPIC, DO_PENTADIAG_INVERSE,       & ! Input !@@ 2p3   6/25/14
          BVPINDEX, BVPSCALEFACTOR, TAYLOR_ORDER, TAYLOR_SMALL,            & ! Input !@@ 2p3/4 6/25/14, 1/7/15
          NLAYERS, NTOTAL, STREAM_VALUE, N_USER_OBSGEOMS, USER_OBSGEOMS,   & ! Inputs     !@@ 2p1
          N_USER_ANGLES, USER_ANGLES, N_USER_RELAZMS, USER_RELAZMS,        & ! Inputs
          FLUX_FACTOR, NBEAMS, BEAM_SZAS, EARTH_RADIUS, HEIGHT_GRID,       & ! Inputs
          DELTAU_INPUT, OMEGA_INPUT, ASYMM_INPUT, D2S_SCALING,             & ! Inputs
          THERMAL_BB_INPUT, LAMBERTIAN_ALBEDO, BRDF_F_0, BRDF_F, UBRDF_F,  & ! Inputs
          EMISSIVITY, SURFBB, SLTERM_ISOTROPIC, SLTERM_F_0,                & ! Inputs  !@@ 2p3 (Add Sleave)
          DO_COLUMN_WFS, DO_SURFACE_WFS, DO_SLEAVE_WFS,                    & ! Inputs  !@@ 2p3 (Add Sleave)
          N_COLUMN_WFS, N_SURFACE_WFS, N_SLEAVE_WFS,                       & ! Inputs  !@@ 2p3 (Add Sleave)
          LSSL_SLTERM_ISOTROPIC, LSSL_SLTERM_F_0,                          & ! Inputs  !@@ 2p3 (New sleave)
          L_DELTAU_INPUT, L_OMEGA_INPUT, L_ASYMM_INPUT, L_D2S_SCALING,     & ! Inputs
          LS_BRDF_F_0, LS_BRDF_F, LS_UBRDF_F, LS_EMISSIVITY,               & ! Inputs
          INTENSITY_TOA, COLUMNWF_TOA, SURFACEWF_TOA,                      & ! Outputs
          INTENSITY_BOA, COLUMNWF_BOA, SURFACEWF_BOA,                      & ! Outputs
          RADLEVEL_UP, RADLEVEL_DN,  N_GEOMETRIES,                         & ! Outputs !@@ 2p2
          COLJACLEVEL_UP, COLJACLEVEL_DN, SURFJACLEVEL_UP, SURFJACLEVEL_DN,& ! Outputs !@@ 2p2
          FLUXES_TOA, COLJACFLUXES_TOA, SURFJACFLUXES_TOA,                 & ! Outputs !@@ 2p3
          FLUXES_BOA, COLJACFLUXES_BOA, SURFJACFLUXES_BOA,                 & ! Outputs !@@ 2p3
          STATUS_INPUTCHECK, C_NMESSAGES, C_MESSAGES, C_ACTIONS,           & ! Exception handling
          STATUS_EXECUTION,  E_MESSAGE, E_TRACE_1, E_TRACE_2 )               ! Exception handling

!  Exception handling

        IF ( STATUS_INPUTCHECK .eq. 1 ) THEN
          write(*,'(a,i4)')'INPUT Check failed from Perturbed Column WF Run # ',T
          write(*,*)' - Number of Messages = ', C_NMESSAGES
          DO k = 1, C_NMESSAGES
            write(*,'(A,I3,A,A)')' - Message # ',K,': ', TRIM(C_MESSAGES(K))
            write(*,'(A,I3,A,A)')' - Action  # ',K,': ', TRIM(C_ACTIONS(K))
          ENDDO
          stop'Test_2S_LPCS program aborted'
        ENDIF
        IF ( STATUS_EXECUTION .eq. 1 ) THEN
          write(*,'(a,i4)')'EXECUTION failed from Perturbed Column WF Run # ',T
          write(*,*)' - Print 1 Message and 2 Traces'
          write(*,'(A)') TRIM(E_MESSAGE)
          write(*,'(A)') TRIM(E_TRACE_1)
          write(*,'(A)') TRIM(E_TRACE_2)
          stop'Test_2S_LPCS program aborted'
        ENDIF

!  Save the last results produced by the OMP thread

        IF (wvn == (TID+1)*numwvn/OMP_NTHREADS) THEN
          Q = 1
          INTENSITY_TOA_SURF_PT(1:N_GEOMETRIES,Q,T,TID+1) = INTENSITY_TOA(1:N_GEOMETRIES)
          INTENSITY_BOA_SURF_PT(1:N_GEOMETRIES,Q,T,TID+1) = INTENSITY_BOA(1:N_GEOMETRIES)

          if ( do_2s_LEVELOUT ) THEN
             RADLEVEL_UP_SURF_PT(1:N_GEOMETRIES,0:NLAYERS,Q,T,TID+1) = RADLEVEL_UP(1:N_GEOMETRIES,0:NLAYERS)
             RADLEVEL_DN_SURF_PT(1:N_GEOMETRIES,0:NLAYERS,Q,T,TID+1) = RADLEVEL_DN(1:N_GEOMETRIES,0:NLAYERS)
          endif

          FLUXES_TOA_SURF_PT(1:NBEAMS,:,Q,T,TID+1) = FLUXES_TOA(1:NBEAMS,:)
          FLUXES_BOA_SURF_PT(1:NBEAMS,:,Q,T,TID+1) = FLUXES_BOA(1:NBEAMS,:)
        END IF

!  End thread loop

      enddo


!  FD calculation 2 : Perturbation on mol abs tau
!  ==============================================

      !write(*,*)
      !write(*,*)'Doing OMP thread ',TID,' Doing FD calculation 2: perturb mol abs tau'

!  Start thread loop

      do thread = 1, nthreads
!      do thread = 3, 3
        t = thread

!  Thread 1: 2S, Solar only   , no   delta-M scaling
!  Thread 2: 2S, Solar only   , with delta-M scaling
!  Thread 3: 2S, Thermal only , no   delta-M scaling
!  Thread 4: 2S, Thermal only , with delta-M scaling
!  Thread 5: 2S, Solar+Thermal, no   delta-M scaling
!  Thread 6: 2S, Solar+Thermal, with delta-M scaling

        call thread_setter &
          ( thread, n_user_obsgeoms, & 
            do_solar_sources, do_thermal_emission, &
            do_d2s_scaling, do_surface_emission,   &
            nbeams, n_user_relazms )
        !write(*,*)'Doing thread # ', thread
        !write(*,'(1x,a,i1,a,i1)') 'Doing OMP thread ',TID,' LIDORT thread ',thread

        column = 1.0d0

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
          !emissivity    =  1.0d0 - albedo_save
        else
          surfbb        = 0.0d0
          emissivity    = 0.0d0
        end if

!  Create optical properties, every thread

        do n = 1, n6
          molsca = molomg(n) * molext(n)
          molabs = molext(n) - molsca

          deltau_input(n) = molsca + molabs*epsfac
          omega_input(n)  = molsca/deltau_input(n)
          asymm_input(n)  = 0.0d0

          m1 = raymoms(2,n) * molsca ; m2 = 0.0d0
          d2s_scaling(n)  = ( m1 + m2 ) / molsca / 5.0d0
        enddo

        parcel = taer / ( height_grid(n6) - height_grid(nlayers) )
        do n = n6 + 1, nlayers
          molsca = molomg(n) * molext(n)
          molabs = molext(n) - molsca

          aerext = Parcel * ( height_grid(n-1) - height_grid(n) )
          aersca = aerext * waer

          totext = molsca + molabs*epsfac + aerext
          totsca = molsca + aersca
          raywt  = molsca / totsca
          aerwt  = aersca / totsca
          omega  = totsca / totext

          deltau_input(n) = totext
          omega_input(n)  = omega
          asymm_input(n)  = gaer * aerwt

          m1 = raymoms(2,n) * molsca
          m2 = 5.0d0 * gaer * gaer * aersca
          d2s_scaling(n)  = ( m1 + m2 ) / totsca / 5.0d0
        enddo

!  Set up surface BRDF

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

          !Linearized BRDF supplement inputs
          DO_KERNEL_FACTOR_WFS      = .FALSE.
          DO_KERNEL_PARAMS_WFS      = .FALSE.

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
          !      TWOSTREAM_LS_BRDFMASTER using the DO_SHADOW_EFFECT flag

          !Call linearized BRDF supplement
          CALL TWOSTREAM_LS_BRDFMASTER &
            ( MAXBEAMS, MAX_USER_ANGLES, MAX_USER_OBSGEOMS,      & ! Dimensions !@@
              MAXSTREAMS_BRDF, MAX_BRDF_KERNELS,                 & ! Dimensions
              MAX_BRDF_PARAMETERS, MAX_SURFACEWFS,               & ! Dimensions
              DO_SOLAR_SOURCES, DO_USER_OBSGEOMS,                & ! Inputs !@@
              LAMBERTIAN_KERNEL_FLAG,                            & ! Inputs
              DO_SHADOW_EFFECT, DO_SURFACE_EMISSION,             & ! Inputs
              NBEAMS, N_USER_ANGLES, N_USER_OBSGEOMS,            & ! Inputs !@@
              BEAM_SZAS, USER_ANGLES, USER_OBSGEOMS,             & ! Inputs !@@
              STREAM_VALUE, NSTREAMS_BRDF,                       & ! Inputs
              N_BRDF_KERNELS, WHICH_BRDF, BRDF_FACTORS,          & ! Inputs
              N_BRDF_PARAMETERS, BRDF_PARAMETERS,                & ! Inputs
              DO_KERNEL_FACTOR_WFS, DO_KERNEL_PARAMS_WFS,        & ! Inputs
              DO_KPARAMS_DERIVS, N_SURFACE_WFS,                  & ! Outputs
              N_KERNEL_FACTOR_WFS, N_KERNEL_PARAMS_WFS,          & ! Outputs
              BRDF_F_0, BRDF_F, UBRDF_F, EMISSIVITY,             & ! Outputs
              LS_BRDF_F_0, LS_BRDF_F, LS_UBRDF_F, LS_EMISSIVITY, & ! Outputs
              STATUS_BRDFSUP, MESSAGE_BRDF, ACTION_BRDF )          ! Outputs

        else
          !Lambertian surface

          LAMBERTIAN_ALBEDO = albedo_save
          emissivity = 1.0d0 - albedo_save
        end if

!  Exception handling on BRDF

        if ( DO_BRDF_SURFACE ) then
           if ( STATUS_BRDFSUP .ne. 0 ) then
              write(*,'(a,i4)')'BRDF Lin supplement Check failed from Perturbed Column WF Run # ',T
              write(*,*)' - Print 1 Message and 1 Action'
              write(*,'(A)') TRIM(MESSAGE_BRDF)
              write(*,'(A)') TRIM(ACTION_BRDF)
              stop'Test_2S_LPCS program aborted'
           endif
        endif

!  Call to Linearized model

        CALL TWOSTREAM_LCS_MASTER &
        ( MAXLAYERS, MAXTOTAL, MAXMESSAGES, MAXBEAMS, MAX_GEOMETRIES,      & ! Dimensions
          MAX_USER_ANGLES, MAX_USER_RELAZMS, MAX_USER_OBSGEOMS,            & ! Dimensions !@@ 2p1
          MAX_ATMOSWFS, MAX_SURFACEWFS, MAX_SLEAVEWFS,                     & ! Dimensions !@@ 2p3 (Add Sleave)
          DO_UPWELLING, DO_DNWELLING, DO_PLANE_PARALLEL, DO_2S_LEVELOUT,   & ! Inputs     !@@ 2p2
          DO_MVOUT_ONLY, DO_ADDITIONAL_MVOUT,                              & ! Inputs     !@@ 2p3
          DO_SOLAR_SOURCES, DO_THERMAL_EMISSION, DO_SURFACE_EMISSION,      & ! Inputs
          DO_D2S_SCALING, DO_BRDF_SURFACE, DO_USER_OBSGEOMS,               & ! Inputs     !@@ 2p1
          DO_SURFACE_LEAVING, DO_SL_ISOTROPIC, DO_PENTADIAG_INVERSE,       & ! Input !@@ 2p3   6/25/14
          BVPINDEX, BVPSCALEFACTOR, TAYLOR_ORDER, TAYLOR_SMALL,            & ! Input !@@ 2p3/4 6/25/14, 1/7/15
          NLAYERS, NTOTAL, STREAM_VALUE, N_USER_OBSGEOMS, USER_OBSGEOMS,   & ! Inputs     !@@ 2p1
          N_USER_ANGLES, USER_ANGLES, N_USER_RELAZMS, USER_RELAZMS,        & ! Inputs
          FLUX_FACTOR, NBEAMS, BEAM_SZAS, EARTH_RADIUS, HEIGHT_GRID,       & ! Inputs
          DELTAU_INPUT, OMEGA_INPUT, ASYMM_INPUT, D2S_SCALING,             & ! Inputs
          THERMAL_BB_INPUT, LAMBERTIAN_ALBEDO, BRDF_F_0, BRDF_F, UBRDF_F,  & ! Inputs
          EMISSIVITY, SURFBB, SLTERM_ISOTROPIC, SLTERM_F_0,                & ! Inputs  !@@ 2p3 (Add Sleave)
          DO_COLUMN_WFS, DO_SURFACE_WFS, DO_SLEAVE_WFS,                    & ! Inputs  !@@ 2p3 (Add Sleave)
          N_COLUMN_WFS, N_SURFACE_WFS, N_SLEAVE_WFS,                       & ! Inputs  !@@ 2p3 (Add Sleave)
          LSSL_SLTERM_ISOTROPIC, LSSL_SLTERM_F_0,                          & ! Inputs  !@@ 2p3 (New sleave)
          L_DELTAU_INPUT, L_OMEGA_INPUT, L_ASYMM_INPUT, L_D2S_SCALING,     & ! Inputs
          LS_BRDF_F_0, LS_BRDF_F, LS_UBRDF_F, LS_EMISSIVITY,               & ! Inputs
          INTENSITY_TOA, COLUMNWF_TOA, SURFACEWF_TOA,                      & ! Outputs
          INTENSITY_BOA, COLUMNWF_BOA, SURFACEWF_BOA,                      & ! Outputs
          RADLEVEL_UP, RADLEVEL_DN,  N_GEOMETRIES,                         & ! Outputs !@@ 2p2
          COLJACLEVEL_UP, COLJACLEVEL_DN, SURFJACLEVEL_UP, SURFJACLEVEL_DN,& ! Outputs !@@ 2p2
          FLUXES_TOA, COLJACFLUXES_TOA, SURFJACFLUXES_TOA,                 & ! Outputs !@@ 2p3
          FLUXES_BOA, COLJACFLUXES_BOA, SURFJACFLUXES_BOA,                 & ! Outputs !@@ 2p3
          STATUS_INPUTCHECK, C_NMESSAGES, C_MESSAGES, C_ACTIONS,           & ! Exception handling
          STATUS_EXECUTION,  E_MESSAGE, E_TRACE_1, E_TRACE_2 )               ! Exception handling

!  Exception handling

        IF ( STATUS_INPUTCHECK .eq. 1 ) THEN
          write(*,'(a,i4)')'INPUT Check failed from Perturbed Column WF Run # ',T
          write(*,*)' - Number of Messages = ', C_NMESSAGES
          DO k = 1, C_NMESSAGES
            write(*,'(A,I3,A,A)')' - Message # ',K,': ', TRIM(C_MESSAGES(K))
            write(*,'(A,I3,A,A)')' - Action  # ',K,': ', TRIM(C_ACTIONS(K))
          ENDDO
          stop'Test_2S_LPCS program aborted'
        ENDIF
        IF ( STATUS_EXECUTION .eq. 1 ) THEN
          write(*,'(a,i4)')'EXECUTION failed from Perturbed Column WF Run # ',T
          write(*,*)' - Print 1 Message and 2 Traces'
          write(*,'(A)') TRIM(E_MESSAGE)
          write(*,'(A)') TRIM(E_TRACE_1)
          write(*,'(A)') TRIM(E_TRACE_2)
          stop'Test_2S_LPCS program aborted'
        ENDIF

!  Save the last results produced by the OMP thread

        IF (wvn == (TID+1)*numwvn/OMP_NTHREADS) THEN
          Q = 1
          INTENSITY_TOA_COL_PT(1:N_GEOMETRIES,Q,T,TID+1) = INTENSITY_TOA(1:N_GEOMETRIES)
          INTENSITY_BOA_COL_PT(1:N_GEOMETRIES,Q,T,TID+1) = INTENSITY_BOA(1:N_GEOMETRIES)

          if ( do_2s_LEVELOUT ) THEN
             RADLEVEL_UP_COL_PT(1:N_GEOMETRIES,0:NLAYERS,Q,T,TID+1) = RADLEVEL_UP(1:N_GEOMETRIES,0:NLAYERS)
             RADLEVEL_DN_COL_PT(1:N_GEOMETRIES,0:NLAYERS,Q,T,TID+1) = RADLEVEL_DN(1:N_GEOMETRIES,0:NLAYERS)
          endif

          FLUXES_TOA_COL_PT(1:NBEAMS,:,Q,T,TID+1) = FLUXES_TOA(1:NBEAMS,:)
          FLUXES_BOA_COL_PT(1:NBEAMS,:,Q,T,TID+1) = FLUXES_BOA(1:NBEAMS,:)
        END IF

!  End thread loop

      enddo


!  FD calculation 3 : Perturbation on aero ext tau
!  ===============================================

      !write(*,*)
      !write(*,*)'Doing OMP thread ',TID,' Doing FD calculation 3: perturb aero ext tau'

!  Start thread loop

      do thread = 1, nthreads
!      do thread = 3, 3
        t = thread

!  Thread 1: 2S, Solar only   , no   delta-M scaling
!  Thread 2: 2S, Solar only   , with delta-M scaling
!  Thread 3: 2S, Thermal only , no   delta-M scaling
!  Thread 4: 2S, Thermal only , with delta-M scaling
!  Thread 5: 2S, Solar+Thermal, no   delta-M scaling
!  Thread 6: 2S, Solar+Thermal, with delta-M scaling

        call thread_setter &
          ( thread, n_user_obsgeoms, & 
            do_solar_sources, do_thermal_emission, &
            do_d2s_scaling, do_surface_emission,   &
            nbeams, n_user_relazms )
        !write(*,*)'Doing thread # ', thread
        !write(*,'(1x,a,i1,a,i1)') 'Doing OMP thread ',TID,' LIDORT thread ',thread

        column = 1.0d0

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
          !emissivity    =  1.0d0 - albedo_save
        else
          surfbb        = 0.0d0
          emissivity    = 0.0d0
        end if

!  Create optical properties, every thread

        do n = 1, n6
          molsca = molomg(n) * molext(n)
          molabs = molext(n) - molsca

          deltau_input(n) = molext(n)
          omega_input(n)  = molomg(n)
          asymm_input(n)  = 0.0d0

          m1     = raymoms(2,n) * molsca ; m2 = 0.0d0
          d2s_scaling(n)  = ( m1 + m2 ) / molsca / 5.0d0
        enddo

        !alternate method
        !taer = taer*epsfac

        parcel = taer / ( height_grid(n6) - height_grid(nlayers) )
        do n = n6 + 1, nlayers
          molsca = molomg(n) * molext(n)
          molabs = molext(n) - molsca

          aerext = Parcel * ( height_grid(n-1) - height_grid(n) )*epsfac
          !aerext = Parcel * ( height_grid(n-1) - height_grid(n) )
          aersca = aerext * waer

          totext = molsca + molabs + aerext
          totsca = molsca + aersca
          raywt  = molsca / totsca
          aerwt  = aersca / totsca
          omega  = totsca / totext

          deltau_input(n) = totext
          omega_input(n)  = omega
          asymm_input(n)  = gaer * aerwt

          m1 = raymoms(2,n) * molsca
          m2 = 5.0d0 * gaer * gaer * aersca
          d2s_scaling(n)  = ( m1 + m2 ) / totsca / 5.0d0
        enddo

!  Set up surface BRDF

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

          !Linearized BRDF supplement inputs
          DO_KERNEL_FACTOR_WFS      = .FALSE.
          DO_KERNEL_PARAMS_WFS      = .FALSE.

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
          !      TWOSTREAM_LS_BRDFMASTER using the DO_SHADOW_EFFECT flag

          !Call linearized BRDF supplement
          CALL TWOSTREAM_LS_BRDFMASTER &
            ( MAXBEAMS, MAX_USER_ANGLES, MAX_USER_OBSGEOMS,      & ! Dimensions !@@
              MAXSTREAMS_BRDF, MAX_BRDF_KERNELS,                 & ! Dimensions
              MAX_BRDF_PARAMETERS, MAX_SURFACEWFS,               & ! Dimensions
              DO_SOLAR_SOURCES, DO_USER_OBSGEOMS,                & ! Inputs !@@
              LAMBERTIAN_KERNEL_FLAG,                            & ! Inputs
              DO_SHADOW_EFFECT, DO_SURFACE_EMISSION,             & ! Inputs
              NBEAMS, N_USER_ANGLES, N_USER_OBSGEOMS,            & ! Inputs !@@
              BEAM_SZAS, USER_ANGLES, USER_OBSGEOMS,             & ! Inputs !@@
              STREAM_VALUE, NSTREAMS_BRDF,                       & ! Inputs
              N_BRDF_KERNELS, WHICH_BRDF, BRDF_FACTORS,          & ! Inputs
              N_BRDF_PARAMETERS, BRDF_PARAMETERS,                & ! Inputs
              DO_KERNEL_FACTOR_WFS, DO_KERNEL_PARAMS_WFS,        & ! Inputs
              DO_KPARAMS_DERIVS, N_SURFACE_WFS,                  & ! Outputs
              N_KERNEL_FACTOR_WFS, N_KERNEL_PARAMS_WFS,          & ! Outputs
              BRDF_F_0, BRDF_F, UBRDF_F, EMISSIVITY,             & ! Outputs
              LS_BRDF_F_0, LS_BRDF_F, LS_UBRDF_F, LS_EMISSIVITY, & ! Outputs
              STATUS_BRDFSUP, MESSAGE_BRDF, ACTION_BRDF )          ! Outputs

        else
          !Lambertian surface

          LAMBERTIAN_ALBEDO = albedo_save
          emissivity = 1.0d0 - albedo_save
        end if

!  Exception handling on BRDF

        if ( DO_BRDF_SURFACE ) then
           if ( STATUS_BRDFSUP .ne. 0 ) then
              write(*,'(a,i4)')'BRDF Lin supplement Check failed from Perturbed Column WF Run # ',T
              write(*,*)' - Print 1 Message and 1 Action'
              write(*,'(A)') TRIM(MESSAGE_BRDF)
              write(*,'(A)') TRIM(ACTION_BRDF)
              stop'Test_2S_LPCS program aborted'
           endif
        endif

!  Call to Linearized model

        CALL TWOSTREAM_LCS_MASTER &
        ( MAXLAYERS, MAXTOTAL, MAXMESSAGES, MAXBEAMS, MAX_GEOMETRIES,      & ! Dimensions
          MAX_USER_ANGLES, MAX_USER_RELAZMS, MAX_USER_OBSGEOMS,            & ! Dimensions !@@ 2p1
          MAX_ATMOSWFS, MAX_SURFACEWFS, MAX_SLEAVEWFS,                     & ! Dimensions !@@ 2p3 (Add Sleave)
          DO_UPWELLING, DO_DNWELLING, DO_PLANE_PARALLEL, DO_2S_LEVELOUT,   & ! Inputs     !@@ 2p2
          DO_MVOUT_ONLY, DO_ADDITIONAL_MVOUT,                              & ! Inputs     !@@ 2p3
          DO_SOLAR_SOURCES, DO_THERMAL_EMISSION, DO_SURFACE_EMISSION,      & ! Inputs
          DO_D2S_SCALING, DO_BRDF_SURFACE, DO_USER_OBSGEOMS,               & ! Inputs     !@@ 2p1
          DO_SURFACE_LEAVING, DO_SL_ISOTROPIC, DO_PENTADIAG_INVERSE,       & ! Input !@@ 2p3   6/25/14
          BVPINDEX, BVPSCALEFACTOR, TAYLOR_ORDER, TAYLOR_SMALL,            & ! Input !@@ 2p3/4 6/25/14, 1/7/15
          NLAYERS, NTOTAL, STREAM_VALUE, N_USER_OBSGEOMS, USER_OBSGEOMS,   & ! Inputs     !@@ 2p1
          N_USER_ANGLES, USER_ANGLES, N_USER_RELAZMS, USER_RELAZMS,        & ! Inputs
          FLUX_FACTOR, NBEAMS, BEAM_SZAS, EARTH_RADIUS, HEIGHT_GRID,       & ! Inputs
          DELTAU_INPUT, OMEGA_INPUT, ASYMM_INPUT, D2S_SCALING,             & ! Inputs
          THERMAL_BB_INPUT, LAMBERTIAN_ALBEDO, BRDF_F_0, BRDF_F, UBRDF_F,  & ! Inputs
          EMISSIVITY, SURFBB, SLTERM_ISOTROPIC, SLTERM_F_0,                & ! Inputs  !@@ 2p3 (Add Sleave)
          DO_COLUMN_WFS, DO_SURFACE_WFS, DO_SLEAVE_WFS,                    & ! Inputs  !@@ 2p3 (Add Sleave)
          N_COLUMN_WFS, N_SURFACE_WFS, N_SLEAVE_WFS,                       & ! Inputs  !@@ 2p3 (Add Sleave)
          LSSL_SLTERM_ISOTROPIC, LSSL_SLTERM_F_0,                          & ! Inputs  !@@ 2p3 (New sleave)
          L_DELTAU_INPUT, L_OMEGA_INPUT, L_ASYMM_INPUT, L_D2S_SCALING,     & ! Inputs
          LS_BRDF_F_0, LS_BRDF_F, LS_UBRDF_F, LS_EMISSIVITY,               & ! Inputs
          INTENSITY_TOA, COLUMNWF_TOA, SURFACEWF_TOA,                      & ! Outputs
          INTENSITY_BOA, COLUMNWF_BOA, SURFACEWF_BOA,                      & ! Outputs
          RADLEVEL_UP, RADLEVEL_DN,  N_GEOMETRIES,                         & ! Outputs !@@ 2p2
          COLJACLEVEL_UP, COLJACLEVEL_DN, SURFJACLEVEL_UP, SURFJACLEVEL_DN,& ! Outputs !@@ 2p2
          FLUXES_TOA, COLJACFLUXES_TOA, SURFJACFLUXES_TOA,                 & ! Outputs !@@ 2p3
          FLUXES_BOA, COLJACFLUXES_BOA, SURFJACFLUXES_BOA,                 & ! Outputs !@@ 2p3
          STATUS_INPUTCHECK, C_NMESSAGES, C_MESSAGES, C_ACTIONS,           & ! Exception handling
          STATUS_EXECUTION,  E_MESSAGE, E_TRACE_1, E_TRACE_2 )               ! Exception handling

!  Exception handling

        IF ( STATUS_INPUTCHECK .eq. 1 ) THEN
          write(*,'(a,i4)')'INPUT Check failed from Perturbed Column WF Run # ',T
          write(*,*)' - Number of Messages = ', C_NMESSAGES
          DO k = 1, C_NMESSAGES
            write(*,'(A,I3,A,A)')' - Message # ',K,': ', TRIM(C_MESSAGES(K))
            write(*,'(A,I3,A,A)')' - Action  # ',K,': ', TRIM(C_ACTIONS(K))
          ENDDO
          stop'Test_2S_LPCS program aborted'
        ENDIF
        IF ( STATUS_EXECUTION .eq. 1 ) THEN
          write(*,'(a,i4)')'EXECUTION failed from Perturbed Column WF Run # ',T
          write(*,*)' - Print 1 Message and 2 Traces'
          write(*,'(A)') TRIM(E_MESSAGE)
          write(*,'(A)') TRIM(E_TRACE_1)
          write(*,'(A)') TRIM(E_TRACE_2)
          stop'Test_2S_LPCS program aborted'
        ENDIF

!  Save the last results produced by the OMP thread

        IF (wvn == (TID+1)*numwvn/OMP_NTHREADS) THEN
          Q = 2
          INTENSITY_TOA_COL_PT(1:N_GEOMETRIES,Q,T,TID+1) = INTENSITY_TOA(1:N_GEOMETRIES)
          INTENSITY_BOA_COL_PT(1:N_GEOMETRIES,Q,T,TID+1) = INTENSITY_BOA(1:N_GEOMETRIES)

          if ( do_2s_LEVELOUT ) THEN
             RADLEVEL_UP_COL_PT(1:N_GEOMETRIES,0:NLAYERS,Q,T,TID+1) = RADLEVEL_UP(1:N_GEOMETRIES,0:NLAYERS)
             RADLEVEL_DN_COL_PT(1:N_GEOMETRIES,0:NLAYERS,Q,T,TID+1) = RADLEVEL_DN(1:N_GEOMETRIES,0:NLAYERS)
          endif

          FLUXES_TOA_COL_PT(1:NBEAMS,:,Q,T,TID+1) = FLUXES_TOA(1:NBEAMS,:)
          FLUXES_BOA_COL_PT(1:NBEAMS,:,Q,T,TID+1) = FLUXES_BOA(1:NBEAMS,:)
        END IF

!  End thread loop

      enddo

!  End "wavenumber" loop

      enddo

!$OMP END DO

      !IF (TID == 0) call cpu_time(omp_e2)

!  End parallel region

!$OMP END PARALLEL

!  Timing tests (OpenMP)
!  =====================

      call cpu_time(omp_e2)

      if (omp_nthreads <= n_core) then
        time_divider = omp_nthreads
      else
        time_divider = n_core
      endif

      write(*,*)
      write(*,'(15x,a)')                    'Timing report'
      write(*,'(4(1x,a))') 'Numwvn','# OMP Threads','Time (sec)'
      write(*,'(4(1x,a))') '------','-------------','----------'
      write(*,'(1x,i5,8x,i1,7x,f6.2)') &
        numwvn, omp_nthreads, (omp_e2 - omp_e1)/real(time_divider)

!  Write Part 2 Results
!  --------------------

      CALL WRITE_2S_LCS_OUTPUT()

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!            END OF MAIN PROFILE/COLUMN/SURFACE WFS TESTS             @
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  Begin file close loop

      do ompt=1,OMP_MAXTHREADS

!  Close output files

        close(10*ompt+36)
        close(10*ompt+37)
        if ( do_2S_LEVELOUT ) then
          close(10*ompt+336)
          close(10*ompt+337)
        endif

!  End file close loop

      enddo

!  End test loop

      enddo

!  Finish

      write(*,*)
      STOP 'successful run'

      CONTAINS

!******************************************************************************
SUBROUTINE WRITE_2S_LPS_OUTPUT()

!  Begin file write loop

      do ompt=1,OMP_MAXTHREADS

!  Write Part 1 results: TOA/BOA Intensity Output
!  ==============================================

!  Start geometry loop
      do v = 1, N_GEOMETRIES_SAVE

!  (1) TOA file
      if (v == 1) then
        write(10*ompt+36,*)
        write(10*ompt+36,'(/T16,a/T16,a/)') &
          'INTENSITY, 1 SURFACE JACOBIAN, 2 PROFILE JACOBIANS', &
          '=================================================='
      end if

      write(10*ompt+36,'(/T16,a,I2,/T16,a/)')' ALL OUTPUT , threads 1-6, GEOMETRY # ',V, &
                                     '========================================='
      write(10*ompt+36,'(T16,6(6x,a16,6x)/)') &
        '   Solar No DM  ','  Solar Yes DM  ', &
        '  Thermal No DM ',' Thermal Yes DM ', &
        ' Solar+Thermal N',' Solar+Thermal Y'

      write(10*ompt+36,366) 'Intensity' , (INTENSITY_TOA_BAS(V,t,ompt),t=1,nthreads)

      write(10*ompt+36,*)
      q=1 !wrt Cox-Munk slope squared parameter
      write(10*ompt+36,367) 'Surface WF', &
                    (BRDF_PAR*surfacewf_TOA_BAS(V,q,t,ompt),&
                     (INTENSITY_TOA_SURF_PT(V,q,t,ompt)-INTENSITY_TOA_BAS(V,t,ompt))/eps,&
                     t=1,nthreads)

      write(10*ompt+36,*)
      do q=1,2 !wrt (1) mol abs tau & (2) aero ext tau
        do n = 1, nlayers
          write(10*ompt+36,368) 'Profile WF # ' ,n, &
                        (profilewf_TOA_BAS(V,n,q,t,ompt),&
                         (INTENSITY_TOA_PROF_PT(V,n,q,t,ompt)-INTENSITY_TOA_BAS(V,t,ompt))/eps,&
                         t=1,nthreads)
        enddo
        write(10*ompt+36,*)
      enddo

!  (2) BOA file
      if (v == 1) then
        write(10*ompt+37,*)
        write(10*ompt+37,'(/T16,a/T16,a/)') &
          'INTENSITY, 1 SURFACE JACOBIAN, 2 PROFILE JACOBIANS', &
          '=================================================='
      end if

      write(10*ompt+37,'(/T16,a,I2,/T16,a/)')' ALL OUTPUT , threads 1-6, GEOMETRY # ',V, &
                                     '========================================='
      write(10*ompt+37,'(T16,6(6x,a16,6x)/)') &
        '   Solar No DM  ','  Solar Yes DM  ', &
        '  Thermal No DM ',' Thermal Yes DM ', &
        ' Solar+Thermal N',' Solar+Thermal Y'

      write(10*ompt+37,366) 'Intensity' , (INTENSITY_BOA_BAS(V,t,ompt),t=1,nthreads)

      write(10*ompt+37,*)
      q=1 !wrt Cox-Munk slope squared parameter
      write(10*ompt+37,367) 'Surface WF', &
                    (BRDF_PAR*surfacewf_BOA_BAS(V,q,t,ompt),&
                     (INTENSITY_BOA_SURF_PT(V,q,t,ompt)-INTENSITY_BOA_BAS(V,t,ompt))/eps,&
                     t=1,nthreads)

      write(10*ompt+37,*)
      do q=1,2 !wrt (1) mol abs tau & (2) aero ext tau
        do n = 1, nlayers
          write(10*ompt+37,368) 'Profile WF # ' ,n, &
                        (profilewf_BOA_BAS(V,n,q,t,ompt),&
                         (INTENSITY_BOA_PROF_PT(V,n,q,t,ompt)-INTENSITY_BOA_BAS(V,t,ompt))/eps,&
                         t=1,nthreads)
        enddo
        write(10*ompt+37,*)
      enddo

!  End geometry loop
      enddo


!  Write Part 1 results: TOA/BOA Flux Output
!  =========================================

!  Start flux loop
      do f = 1, 2 !actinic vs regular
!  Start beam loop
      do v = 1, nbeams

!  (1) TOA file
      write(10*ompt+36,*)
      if ((f == 1) .and. (v == 1)) then
        write(10*ompt+36,'(/T16,a/T16,a/)') &
          'ACTINIC FLUX, 1 SURFACE JACOBIAN, 2 PROFILE JACOBIANS', &
          '====================================================='
      else if ((f == 2) .and. (v == 1)) then
        write(10*ompt+36,'(/T16,a/T16,a/)') &
          'REGULAR FLUX, 1 SURFACE JACOBIAN, 2 PROFILE JACOBIANS', &
          '====================================================='
      end if

      write(10*ompt+36,'(/T16,a,I2,/T16,a/)')' ALL OUTPUT , threads 1-2, BEAM # ',V, &
                                     '====================================='
      write(10*ompt+36,'(T16,2(6x,a16,6x)/)') &
        '   Solar No DM  ','  Solar Yes DM  '

      write(10*ompt+36,366) 'Flux' , (FLUXES_TOA_BAS(V,f,t,ompt),t=1,2)

      write(10*ompt+36,*)
      q=1 !wrt Cox-Munk slope squared parameter
      write(10*ompt+36,367) 'Surface WF', &
                    ( BRDF_PAR*SURFJACFLUXES_TOA_BAS(V,f,q,t,ompt),&
                      (FLUXES_TOA_SURF_PT(V,f,q,t,ompt)-FLUXES_TOA_BAS(V,f,t,ompt))/eps,t=1,2 )

      write(10*ompt+36,*)
      do q=1,2 !wrt (1) mol abs tau & (2) aero ext tau
        do n = 1, nlayers
          write(10*ompt+36,368) 'Profile WF # ' ,n, &
                        ( PROFJACFLUXES_TOA_BAS(V,f,n,q,t,ompt),&
                          (FLUXES_TOA_PROF_PT(V,f,n,q,t,ompt)-FLUXES_TOA_BAS(V,f,t,ompt))/eps,t=1,2 )
        enddo
        write(10*ompt+36,*)
      enddo

!  (2) BOA file
      write(10*ompt+37,*)
      if ((f == 1) .and. (v == 1)) then
        write(10*ompt+37,'(/T16,a/T16,a/)') &
          'ACTINIC FLUX, 1 SURFACE JACOBIAN, 2 PROFILE JACOBIANS', &
          '====================================================='
      else if ((f == 2) .and. (v == 1)) then
        write(10*ompt+37,'(/T16,a/T16,a/)') &
          'REGULAR FLUX, 1 SURFACE JACOBIAN, 2 PROFILE JACOBIANS', &
          '====================================================='
      end if

      write(10*ompt+37,'(/T16,a,I2,/T16,a/)')' ALL OUTPUT , threads 1-2, BEAM # ',V, &
                                     '====================================='
      write(10*ompt+37,'(T16,2(6x,a16,6x)/)') &
        '   Solar No DM  ','  Solar Yes DM  '

      write(10*ompt+37,366) 'Flux' , (FLUXES_BOA_BAS(V,f,t,ompt),t=1,2)

      write(10*ompt+37,*)
      q=1 !wrt Cox-Munk slope squared parameter
      write(10*ompt+37,367) 'Surface WF', &
                    ( BRDF_PAR*SURFJACFLUXES_BOA_BAS(V,f,q,t,ompt),&
                      (FLUXES_BOA_SURF_PT(V,f,q,t,ompt)-FLUXES_BOA_BAS(V,f,t,ompt))/eps,t=1,2 )

      write(10*ompt+37,*)
      do q=1,2 !wrt (1) mol abs tau & (2) aero ext tau
        do n = 1, nlayers
          write(10*ompt+37,368) 'Profile WF # ' ,n, &
                        ( PROFJACFLUXES_BOA_BAS(V,f,n,q,t,ompt),&
                          (FLUXES_BOA_PROF_PT(V,f,n,q,t,ompt)-FLUXES_BOA_BAS(V,f,t,ompt))/eps,t=1,2 )
        enddo
        write(10*ompt+37,*)
      enddo

!  End beam loop
      enddo
!  End flux loop
      enddo

366   format(a,T16,6(7x,e14.6,7x))
367   format(a,T16,6(2e14.6))
368   format(a,i2,T16,6(2e14.6))

!  Write Part 1 results: LEVOUT Intensity Output
!  =============================================

      if ( do_2S_LEVELOUT ) then

!  Start geometry loop
      do v = 1, N_GEOMETRIES_SAVE

!  (1) TOA file
      if (v == 1) then
        write(10*ompt+336,*)
        write(10*ompt+336,'(/T16,a/T16,a/)') &
          'INTENSITY, 1 SURFACE JACOBIAN, 2 PROFILE JACOBIANS', &
          '=================================================='
      end if

      write(10*ompt+336,'(/T33,a,I2,/T33,a/)')' ALL OUTPUT , threads 1-6, GEOMETRY # ',V, &
                                     '========================================='
      write(10*ompt+336,'(T33,6(6x,a16,6x)/)') &
        '   Solar No DM  ','  Solar Yes DM  ', &
        '  Thermal No DM ',' Thermal Yes DM ', &
        ' Solar+Thermal N',' Solar+Thermal Y'

      do n = 0, nlayers
         write(10*ompt+336,3366)'Intensity @ level # ',DBLE(N), (RADLEVEL_UP_BAS(V,N,T,ompt),T=1,NTHREADS)
      enddo

      write(10*ompt+336,*)
      q=1 !wrt Cox-Munk slope squared parameter
      do n = 0, nlayers
         write(10*ompt+336,3367) 'Surface WF @ level #', DBLE(N), &
                    (BRDF_PAR*SURFJACLEVEL_UP_BAS(V,n,q,t,ompt),&
                     (RADLEVEL_UP_SURF_PT(V,n,q,t,ompt)-RADLEVEL_UP_BAS(V,n,t,ompt))/eps,&
                     t=1,nthreads)
      enddo

      write(10*ompt+336,*)
      do n = 0, nlayers
         do q=1,2 !wrt (1) mol abs tau & (2) aero ext tau
            do k = 1, nlayers
               write(10*ompt+336,3368) 'Profile WF # ',q,' @ level #',k,DBLE(N), &
                        (PROFJACLEVEL_UP_BAS(V,n,k,q,t,ompt),&
                         (RADLEVEL_UP_PROF_PT(V,n,k,q,t,ompt)-RADLEVEL_UP_BAS(V,n,t,ompt))/eps,&
                         t=1,nthreads)
            enddo
            write(10*ompt+336,*)
         enddo
      enddo

!  (2) BOA file
      if (v == 1) then
        write(10*ompt+337,*)
        write(10*ompt+337,'(/T16,a/T16,a/)') &
          'INTENSITY, 1 SURFACE JACOBIAN, 2 PROFILE JACOBIANS', &
          '=================================================='
      end if

      write(10*ompt+337,'(/T33,a,I2,/T33,a/)')' ALL OUTPUT , threads 1-6, GEOMETRY # ',V, &
                                     '========================================='
      write(10*ompt+337,'(T33,6(6x,a16,6x)/)') &
        '   Solar No DM  ','  Solar Yes DM  ', &
        '  Thermal No DM ',' Thermal Yes DM ', &
        ' Solar+Thermal N',' Solar+Thermal Y'

      do n = 0, nlayers
         write(10*ompt+337,3366)'Intensity @ level #',DBLE(N), (RADLEVEL_DN_BAS(V,N,T,ompt),T=1,NTHREADS)
      enddo

      write(10*ompt+337,*)
      q=1 !wrt Cox-Munk slope squared parameter
      do n = 0, nlayers
         write(10*ompt+337,3367) 'Surface WF @ level #', DBLE(N), &
                    (BRDF_PAR*SURFJACLEVEL_DN_BAS(V,n,q,t,ompt),&
                     (RADLEVEL_DN_SURF_PT(V,n,q,t,ompt)-RADLEVEL_DN_BAS(V,n,t,ompt))/eps,&
                     t=1,nthreads)
      enddo

      write(10*ompt+337,*)
      do n = 0, nlayers
         do q=1,2 !wrt (1) mol abs tau & (2) aero ext tau
            do k = 1, nlayers
               write(10*ompt+337,3368) 'Profile WF # ',q,' @ level #',k,DBLE(N), &
                        (PROFJACLEVEL_DN_BAS(V,n,k,q,t,ompt),&
                         (RADLEVEL_DN_PROF_PT(V,n,k,q,t,ompt)-RADLEVEL_DN_BAS(V,n,t,ompt))/eps,&
                         t=1,nthreads)
            enddo
            write(10*ompt+337,*)
         enddo
      enddo

!  End geometry loop
      enddo

3366  format(a,f4.1,T33,6(7x,e14.6,7x))
3367  format(a,f4.1,T33,6(2e14.6))
3368  format(2(a,i2),1x,f4.1,T33,6(2e14.6))

!  End do_2S_LEVELOUT if block
      end if

!  End file write loop

      enddo

END SUBROUTINE WRITE_2S_LPS_OUTPUT

!******************************************************************************
SUBROUTINE WRITE_2S_LCS_OUTPUT()

!  Begin file write loop

      do ompt=1,OMP_MAXTHREADS

!  Write Part 2 results: TOA/BOA Intensity Output
!  ==============================================

!  Start geometry loop
      do v = 1, N_GEOMETRIES_SAVE

!  (1) TOA file
      if (v == 1) then
        write(10*ompt+36,*)
        write(10*ompt+36,*) '************************************************************************'
        write(10*ompt+36,'(/T16,a/T16,a/)') &
          'INTENSITY, 1 SURFACE JACOBIAN, 2 COLUMN JACOBIANS', &
          '=================================================='
      end if

      write(10*ompt+36,'(/T16,a,I2,/T16,a/)')' ALL OUTPUT , threads 1-6, GEOMETRY # ',V, &
                                     '========================================='
      write(10*ompt+36,'(T16,6(6x,a16,6x)/)') &
        '   Solar No DM  ','  Solar Yes DM  ', &
        '  Thermal No DM ',' Thermal Yes DM ', &
        ' Solar+Thermal N',' Solar+Thermal Y'

      write(10*ompt+36,366) 'Intensity' , (INTENSITY_TOA_BAS(V,t,ompt),t=1,nthreads)

      write(10*ompt+36,*)
      q=1 !wrt Cox-Munk slope squared parameter
      write(10*ompt+36,367) 'Surface WF', &
                    (BRDF_PAR*surfacewf_TOA_BAS(V,q,t,ompt),&
                     (INTENSITY_TOA_SURF_PT(V,q,t,ompt)-INTENSITY_TOA_BAS(V,t,ompt))/eps,&
                     t=1,nthreads)

      write(10*ompt+36,*)
      do q=1,2 !wrt (1) mol abs tau & (2) aero ext tau
        write(10*ompt+36,368) 'Column WF # ' ,q, &
                      (columnwf_TOA_BAS(V,q,t,ompt),&
                       (INTENSITY_TOA_COL_PT(V,q,t,ompt)-INTENSITY_TOA_BAS(V,t,ompt))/eps,&
                       t=1,nthreads)
      enddo

!  (2) BOA file
      if (v == 1) then
        write(10*ompt+37,*)
        write(10*ompt+37,*) '************************************************************************'
        write(10*ompt+37,'(/T16,a/T16,a/)') &
          'INTENSITY, 1 SURFACE JACOBIAN, 2 COLUMN JACOBIANS', &
          '=================================================='
      end if

      write(10*ompt+37,'(/T16,a,I2,/T16,a/)')' ALL OUTPUT , threads 1-6, GEOMETRY # ',V, &
                                     '========================================='
      write(10*ompt+37,'(T16,6(6x,a16,6x)/)') &
        '   Solar No DM  ','  Solar Yes DM  ', &
        '  Thermal No DM ',' Thermal Yes DM ', &
        ' Solar+Thermal N',' Solar+Thermal Y'

      write(10*ompt+37,366) 'Intensity' , (INTENSITY_BOA_BAS(V,t,ompt),t=1,nthreads)

      write(10*ompt+37,*)
      q=1 !wrt Cox-Munk slope squared parameter
      write(10*ompt+37,367) 'Surface WF', &
                    (BRDF_PAR*surfacewf_BOA_BAS(V,q,t,ompt),&
                     (INTENSITY_BOA_SURF_PT(V,q,t,ompt)-INTENSITY_BOA_BAS(V,t,ompt))/eps,&
                     t=1,nthreads)

      write(10*ompt+37,*)
      do q=1,2 !wrt (1) mol abs tau & (2) aero ext tau
        write(10*ompt+37,368) 'Column WF # ' ,q, &
                      (columnwf_BOA_BAS(V,q,t,ompt),&
                       (INTENSITY_BOA_COL_PT(V,q,t,ompt)-INTENSITY_BOA_BAS(V,t,ompt))/eps,&
                       t=1,nthreads)
      enddo

!  End geometry loop
      enddo

!  Write Part 2 results: TOA/BOA Flux Output
!  =========================================

!  Start flux loop
      do f = 1, 2 !actinic vs regular
!  Start beam loop
      do v = 1, nbeams

!  (1) TOA file
      write(10*ompt+36,*)
      if ((f == 1) .and. (v == 1)) then
        write(10*ompt+36,'(/T16,a/T16,a/)') &
          'ACTINIC FLUX, 1 SURFACE JACOBIAN, 2 COLUMN JACOBIANS', &
          '===================================================='
      else if ((f == 2) .and. (v == 1)) then
        write(10*ompt+36,'(/T16,a/T16,a/)') &
          'REGULAR FLUX, 1 SURFACE JACOBIAN, 2 COLUMN JACOBIANS', &
          '===================================================='
      end if

      write(10*ompt+36,'(/T16,a,I2,/T16,a/)')' ALL OUTPUT , threads 1-2, BEAM # ',V, &
                                     '====================================='
      write(10*ompt+36,'(T16,2(6x,a16,6x)/)') &
        '   Solar No DM  ','  Solar Yes DM  '

      write(10*ompt+36,366) 'Flux' , (FLUXES_TOA_BAS(V,f,t,ompt),t=1,2)

      write(10*ompt+36,*)
      q=1 !wrt Cox-Munk slope squared parameter
      write(10*ompt+36,367) 'Surface WF', &
                    ( BRDF_PAR*SURFJACFLUXES_TOA_BAS(V,f,q,t,ompt),&
                      (FLUXES_TOA_SURF_PT(V,f,q,t,ompt)-FLUXES_TOA_BAS(V,f,t,ompt))/eps,t=1,2 )

      write(10*ompt+36,*)
      do q=1,2 !wrt (1) mol abs tau & (2) aero ext tau
          write(10*ompt+36,368) 'Column WF # ' ,q,&
                        ( COLJACFLUXES_TOA_BAS(V,f,q,t,ompt),&
                          (FLUXES_TOA_COL_PT(V,f,q,t,ompt)-FLUXES_TOA_BAS(V,f,t,ompt))/eps,t=1,2 )
      enddo

!  (2) BOA file
      write(10*ompt+37,*)
      if ((f == 1) .and. (v == 1)) then
        write(10*ompt+37,'(/T16,a/T16,a/)') &
          'ACTINIC FLUX, 1 SURFACE JACOBIAN, 2 COLUMN JACOBIANS', &
          '===================================================='
      else if ((f == 2) .and. (v == 1)) then
        write(10*ompt+37,'(/T16,a/T16,a/)') &
          'REGULAR FLUX, 1 SURFACE JACOBIAN, 2 COLUMN JACOBIANS', &
          '===================================================='
      end if

      write(10*ompt+37,'(/T16,a,I2,/T16,a/)')' ALL OUTPUT , threads 1-2, BEAM # ',V, &
                                     '====================================='
      write(10*ompt+37,'(T16,2(6x,a16,6x)/)') &
        '   Solar No DM  ','  Solar Yes DM  '

      write(10*ompt+37,366) 'Flux' , (FLUXES_BOA_BAS(V,f,t,ompt),t=1,2)

      write(10*ompt+37,*)
      q=1 !wrt Cox-Munk slope squared parameter
      write(10*ompt+37,367) 'Surface WF', &
                    ( BRDF_PAR*SURFJACFLUXES_BOA_BAS(V,f,q,t,ompt),&
                      (FLUXES_BOA_SURF_PT(V,f,q,t,ompt)-FLUXES_BOA_BAS(V,f,t,ompt))/eps,t=1,2 )

      write(10*ompt+37,*)
      do q=1,2 !wrt (1) mol abs tau & (2) aero ext tau
        write(10*ompt+37,368) 'Column WF # ',q, &
                      ( COLJACFLUXES_BOA_BAS(V,f,q,t,ompt),&
                        (FLUXES_BOA_COL_PT(V,f,q,t,ompt)-FLUXES_BOA_BAS(V,f,t,ompt))/eps,t=1,2 )
      enddo

!  End beam loop
      enddo
!  End flux loop
      enddo

366   format(a,T16,6(7x,e14.6,7x))
367   format(a,T16,6(2e14.6))
368   format(a,i2,T16,6(2e14.6))


!  Write Part 2 results: LEVOUT Intensity Output
!  =============================================

      if ( do_2S_LEVELOUT ) then

!  Start geometry loop
      do v = 1, N_GEOMETRIES_SAVE

!  (1) TOA file
      if (v == 1) then
        write(10*ompt+336,*)
        write(10*ompt+336,*) '************************************************************************'
        write(10*ompt+336,'(/T16,a/T16,a/)') &
          'INTENSITY, 1 SURFACE JACOBIAN, 2 COLUMN JACOBIANS', &
          '=================================================='
      end if

      write(10*ompt+336,'(/T33,a,I2,/T33,a/)')' ALL OUTPUT , threads 1-6, GEOMETRY # ',V, &
                                     '========================================='
      write(10*ompt+336,'(T33,6(6x,a16,6x)/)') &
        '   Solar No DM  ','  Solar Yes DM  ', &
        '  Thermal No DM ',' Thermal Yes DM ', &
        ' Solar+Thermal N',' Solar+Thermal Y'

      do n = 0, nlayers
         write(10*ompt+336,3366)'Level Intensity @ level #',DBLE(N), (RADLEVEL_UP_BAS(V,N,T,ompt),T=1,NTHREADS)
      enddo

      write(10*ompt+336,*)
      q=1 !wrt Cox-Munk slope squared parameter
      do n = 0, nlayers
         write(10*ompt+336,3367) 'Surface WF @ level #', DBLE(N), &
                    (BRDF_PAR*SURFJACLEVEL_UP_BAS(V,n,q,t,ompt),&
                     (RADLEVEL_UP_SURF_PT(V,n,q,t,ompt)-RADLEVEL_UP_BAS(V,n,t,ompt))/eps,&
                     t=1,nthreads)
      enddo

      write(10*ompt+336,*)
      do n = 0, nlayers
         do q=1,2 !wrt (1) mol abs tau & (2) aero ext tau
            write(10*ompt+336,3368) 'Column WF ',q,' @ level #', DBLE(N), &
                        (COLJACLEVEL_UP_BAS(V,n,q,t,ompt),&
                         (RADLEVEL_UP_COL_PT(V,n,q,t,ompt)-RADLEVEL_UP_BAS(V,n,t,ompt))/eps,&
                         t=1,nthreads)
         enddo
      enddo

!  (2) BOA file
      if (v == 1) then
        write(10*ompt+337,*)
        write(10*ompt+337,*) '************************************************************************'
        write(10*ompt+337,'(/T16,a/T16,a/)') &
          'INTENSITY, 1 SURFACE JACOBIAN, 2 COLUMN JACOBIANS', &
          '=================================================='
      end if

      write(10*ompt+337,'(/T33,a,I2,/T33,a/)')' ALL OUTPUT , threads 1-6, GEOMETRY # ',V, &
                                     '========================================='
      write(10*ompt+337,'(T33,6(6x,a16,6x)/)') &
        '   Solar No DM  ','  Solar Yes DM  ', &
        '  Thermal No DM ',' Thermal Yes DM ', &
        ' Solar+Thermal N',' Solar+Thermal Y'

      do n = 0, nlayers
         write(10*ompt+337,3366)'Intensity @ level #',DBLE(N), (RADLEVEL_DN_BAS(V,N,T,ompt),T=1,NTHREADS)
      enddo

      write(10*ompt+337,*)
      q=1 !wrt Cox-Munk slope squared parameter
      do n = 0, nlayers
         write(10*ompt+337,3367) 'Surface WF @ level #', DBLE(N), &
                    (BRDF_PAR*SURFJACLEVEL_DN_BAS(V,n,q,t,ompt),&
                     (RADLEVEL_DN_SURF_PT(V,n,q,t,ompt)-RADLEVEL_DN_BAS(V,n,t,ompt))/eps,&
                     t=1,nthreads)
      enddo

      write(10*ompt+337,*)
      do n = 0, nlayers
         do q=1,2 !wrt (1) mol abs tau & (2) aero ext tau
            write(10*ompt+337,3368) 'Column WF ',q,' @ level #' ,DBLE(N), &
                        (COLJACLEVEL_DN_BAS(V,n,q,t,ompt),&
                         (RADLEVEL_DN_COL_PT(V,n,q,t,ompt)-RADLEVEL_DN_BAS(V,n,t,ompt))/eps,&
                         t=1,nthreads)
         enddo
      enddo

!  End geometry loop
      enddo

3366  format(a,f4.1,T33,6(7x,e14.6,7x))
3367  format(a,f4.1,T33,6(2e14.6))
3368  format(a,i2,a,1x,f4.1,T33,6(2e14.6))

!  End do_2S_LEVELOUT if block
      end if

!  End file write loop

      enddo

END SUBROUTINE WRITE_2S_LCS_OUTPUT

!******************************************************************************

END PROGRAM test_2s_LPCS_2p4_OMP

!******************************************************************************
!******************************************************************************
subroutine thread_setter &
   ( thread, n_user_obsgeoms, & 
     do_solar_sources, do_thermal_emission, &
     do_d2s_scaling, do_surface_emission,   &
     nbeams, n_user_relazms )

  implicit none

  integer, intent(in)    :: thread, n_user_obsgeoms
  logical, intent(inout) :: do_solar_sources, do_thermal_emission
  logical, intent(inout) :: do_d2s_scaling, do_surface_emission
  integer, intent(inout) :: nbeams, n_user_relazms

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
          nbeams = n_user_obsgeoms;  n_user_relazms = n_user_obsgeoms ! Must set these 2, otherwise doesn't work
        else if ( thread .eq. 2 ) then
          DO_SOLAR_SOURCES    = .true.
          DO_D2S_SCALING      = .true.
          DO_THERMAL_EMISSION = .false.
          DO_SURFACE_EMISSION = .false.
          nbeams = n_user_obsgeoms;  n_user_relazms = n_user_obsgeoms ! Must set these 2, otherwise doesn't work
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

return
end subroutine thread_setter

!**********************************************************************
!**********************************************************************
