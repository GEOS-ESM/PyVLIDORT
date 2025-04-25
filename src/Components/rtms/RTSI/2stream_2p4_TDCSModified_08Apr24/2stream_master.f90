! ###########################################################
! #                                                         #
! #             THE TWOSTREAM LIDORT MODEL                  #
! #                                                         #
! #      (LInearized Discrete Ordinate Radiative Transfer)  #
! #       --         -        -        -         -          #
! #                                                         #
! ###########################################################

! ###########################################################
! #                                                         #
! #  Authors :      Robert. J. D. Spurr (1)                 #
! #                 Vijay Natraj        (2)                 #
! #                                                         #
! #  Address (1) :     RT Solutions, Inc.                   #
! #                    9 Channing Street                    #
! #                    Cambridge, MA 02138, USA             #
! #  Tel:             (617) 492 1183                        #
! #  Email :           rtsolutions@verizon.net              #
! #                                                         #
! #  Address (2) :     CalTech                              #
! #                    Department of Planetary Sciences     #
! #                    1200 East California Boulevard       #
! #                    Pasadena, CA 91125                   #
! #  Tel:             (626) 395 6962                        #
! #  Email :           vijay@gps.caltech.edu                #
! #                                                         #
! #  Version 1.0-1.3 :                                      #
! #     Mark 1: October  2010                               #
! #     Mark 2: May      2011, with BRDFs                   #
! #     Mark 3: October  2011, with Thermal sources         #
! #                                                         #
! #  Version 2.0-2.1 :                                      #
! #     Mark 4: November 2012, LCS/LPS Split, Fixed Arrays  #
! #     Mark 5: December 2012, Observation Geometry option  #
! #                                                         #
! #  Version 2.2-2.3 :                                      #
! #     Mark 6: July     2013, Level outputs + control      #
! #     Mark 7: December 2013, Flux outputs  + control      #
! #     Mark 8: January  2014, Surface Leaving + control    #
! #     Mark 9: June     2014, Inverse Pentadiagonal        #
! #                                                         #
! #  Version 2.4 :                                          #
! #     Mark 10: August  2014, Green's function Regular     #
! #     Mark 11: January 2015, Green's function Linearized  #
! #                            Taylor, dethreaded, OpenMP   #
! #                                                         #
! ###########################################################

! #############################################################
! #                                                           #
! #   This Version of LIDORT-2STREAM comes with a GNU-style   #
! #   license. Please read the license carefully.             #
! #                                                           #
! #############################################################

! ###############################################################
! #                                                             #
! # Subroutines in this Module                                  #
! #                                                             #
! #            TWOSTREAM_MASTER (top-level master)              #
! #            TWOSTREAM_FOURIER_MASTER                         #
! #                                                             #
! ###############################################################

module twostream_master_m

Use twostream_inputs_m
Use twostream_writemodules_m
Use twostream_miscsetups_m
Use twostream_solutions_m
Use twostream_bvproblem_m ! tridiagonal routines include 4/22/17, 5/24/17.
Use twostream_intensity_m
Use twostream_thermalsup_m

PUBLIC

contains

!  Rob Fix 5/24/17 for 2p4a, reintroduce classical PI methods
!  Rob Fix 5/24/17 for 2p4a, Develop handles for Conservative scattering
!  Rob Fix 5/24/17 for 2p4a, Flux output at all levels

!  5/8/24. Introduce flag DO_HEMISPHER_MEAN for thermal regime
!    -- Phasefunction is 1 + G for forward  scattering hemisphere
!    -- Phasefunction is 1 - G for backward scattering hemisphere
!    -- Mubar must be 0.5, No Deltam-scaling. Not allowed with solar scattering

!  5/9/24. Introduce flag DO_EDDINGTON_FLUX for solar regime, FLUXES-ONLY
!    -- Mubar must be 0.5, No Deltam-scaling. Not allowed with thermal sources

SUBROUTINE TWOSTREAM_MASTER ( DO_DEBUG_INPUT, &
          MAXLAYERS, MAXTOTAL, MAXMESSAGES, MAXBEAMS, MAX_GEOMETRIES,     & ! Dimensions
          MAX_USER_RELAZMS, MAX_USER_STREAMS, MAX_USER_OBSGEOMS,          & ! Dimensions !@@ 2p1
          DO_UPWELLING, DO_DNWELLING, DO_PLANE_PARALLEL, DO_2S_LEVELOUT,  & ! Inputs     !@@ 2p2
          DO_MVOUT_ONLY, DO_ADDITIONAL_MVOUT, DO_TRIDIAG, DO_CLASSICAL_PI,& ! Inputs     !@@ 2p3
          DO_SOLAR_SOURCES, DO_THERMAL_EMISSION, DO_SURFACE_EMISSION,     & ! Inputs
          DO_HEMISPHER_MEAN, DO_EDDINGTON_FLUX,                           & ! Inputs
          DO_D2S_SCALING, DO_BRDF_SURFACE, DO_USER_OBSGEOMS, CSTHRESH,    & ! Inputs     !@@ 2p1
          DO_SURFACE_LEAVING, DO_SL_ISOTROPIC, DO_PENTADIAG_INVERSE,      & ! Input !@@ 2p3 6/25/14
          BVPINDEX, BVPSCALEFACTOR, TAYLOR_ORDER, TAYLOR_SMALL, TCUTOFF,  & ! Input !@@ 2p3 6/25/14, 8/15/14
          NLAYERS, NTOTAL, STREAM_VALUE, N_USER_OBSGEOMS, USER_OBSGEOMS,  & ! Inputs     !@@ 2p1
          N_USER_STREAMS, USER_ANGLES, N_USER_RELAZMS, USER_RELAZMS,      & ! Inputs
          FLUX_FACTOR, NBEAMS, BEAM_SZAS, EARTH_RADIUS, HEIGHT_GRID,      & ! Inputs
          DELTAU_INPUT, OMEGA_INPUT, ASYMM_INPUT, D2S_SCALING,            & ! Inputs
          THERMAL_BB_INPUT, LAMBERTIAN_ALBEDO, BRDF_F_0, BRDF_F, UBRDF_F, & ! Inputs
          EMISSIVITY, SURFBB, SLTERM_ISOTROPIC, SLTERM_F_0,               & ! Inputs  !@@ 2p3 (Sleave)
          INTENSITY_TOA, INTENSITY_BOA, FLUXES_TOA, FLUXES_BOA,           & ! Outputs !@@ 2p3 (Fluxes)
          RADLEVEL_UP, RADLEVEL_DN, FLUXLEVEL_UP, FLUXLEVEL_DN, N_GEOMETRIES, & ! Outputs !@@ 2p2
          STATUS_INPUTCHECK, C_NMESSAGES, C_MESSAGES, C_ACTIONS,          & ! Exception handling
          STATUS_EXECUTION,  E_MESSAGE, E_TRACE_1, E_TRACE_2 )              ! Exception handling

      implicit none

!  precision and parameters

      INTEGER      , PARAMETER :: dp   = KIND( 1.0D0 )
      REAL(kind=dp), parameter :: zero = 0.0_dp, one = 1.0_dp

!  Notes 21 december 2012. Observational Geometry Inputs. Marked with !@@ 2p1

!     Observation-Geometry New dimensioning.    MAX_USER_OBSGEOMS
!     Observation-Geometry input control.       DO_USER_OBSGEOMS
!     Observation-Geometry input control.       N_USER_OBSGEOMS
!     User-defined Observation Geometry angles. USER_OBSGEOMS

!  Notes 17 July 2013, Optional output at all levels. Marked with !@@ 2p2
!      New flag for input : DO_2S_LEVELOUT

!  Notes 05 November 2013. Flux output options. Two New Flags
!   DO_MVOUT_ONLY
!   DO_ADDITIONAL_MVOUT

!  Notes 23 January 2014. surface Leaving options. Two New Flags
!   DO_SURFACE_LEAVING
!   DO_SL_ISOTROPIC

!  Notes 25 June 2014. BVProblem control
!   * PentaDiagonal Inverse flag (BVP solved from bottom to top). Only for BVPIndex = 1
!   * BVP Index : 0 = LAPACK, 1 = Penta # 1 (original), 2 = Penta # 2 (new, 2012 Kanal paper)
!   * BVP Scale Factor. Debug only. Set this to 1.0 on input

!  Version 2.4 Notes. 15 August 2014. Greens Function implementation.
!  ------------------------------------------------------------------

!   * Greens function particular integral and postprocessing
!   * Use of PPSTREAM maks to reduce coding for observational geometry
!   * Dethreading (removal of MAXTHREADS dimension)
!   * Use of Taylor series routines for limiting cases

!  4/12/24. Introduce CSTHRESH = Threshold for using Conservative scattering case

!  Subroutine input arguments
!  --------------------------

!  5/8/24. Introduce flag DO_HEMISPHER_MEAN for thermal regime

      Logical, intent(in) :: DO_HEMISPHER_MEAN

!  5/9/24. Introduce flag DO_EDDINGTON_FLUX for solar flux regime

      Logical, intent(in) :: DO_EDDINGTON_FLUX

!  debug flag

      Logical, intent(in) :: do_debug_input

!  Dimensions :
!      MAXTOTAL       = 2 * MAXLAYERS
!      MAX_GEOMETRIES = MAXBEAMS * MAX_USER_STREAMS * MAX_USER_RELAZMS
!      !@@ MAX_USER_OBSGEOMS >/= MAXBEAMS

!  @@ Rob Spurr, 15 August 2014, Version 2.4, MAXTHREADS dimension removed

      INTEGER, INTENT(IN)        :: MAXMESSAGES, MAXLAYERS, MAXTOTAL
      INTEGER, INTENT(IN)        :: MAXBEAMS, MAX_GEOMETRIES
      INTEGER, INTENT(IN)        :: MAX_USER_STREAMS, MAX_USER_RELAZMS, MAX_USER_OBSGEOMS

!  Directional Flags

      LOGICAL, INTENT(IN)        :: DO_UPWELLING, DO_DNWELLING

!  MS-only flag (Mode of operation). NOT REQUIRED
!    IF set, only calculating  MS field
!      LOGICAL, INTENT(IN)        :: DO_MSMODE_2STREAM

!  Plane parallel flag

      LOGICAL, INTENT(IN)        :: DO_PLANE_PARALLEL

!  @@ Rob Spurr, 17 July 2013, Version 2.2, Levelout flag

      LOGICAL, INTENT(IN)        :: DO_2S_LEVELOUT     ! @@ 2p2

!  @@ Rob Spurr, 05 November 2013, Version 2.3, Flux option flags

      LOGICAL, INTENT(IN)        :: DO_MVOUT_ONLY       ! @@ 2p3
      LOGICAL, INTENT(IN)        :: DO_ADDITIONAL_MVOUT ! @@ 2p3

!  New Tridiagonal flag

      LOGICAL, INTENT(IN)        :: DO_TRIDIAG

!  Classical PI flag
!  Rob Fix 5/24/17 for 2p4a, reintroduce classical PI methods

      LOGICAL, INTENT(IN)        :: DO_CLASSICAL_PI

!  ** New **. October 2011, Sources control, including thermal

      LOGICAL, INTENT(IN)        :: DO_SOLAR_SOURCES
      LOGICAL, INTENT(IN)        :: DO_THERMAL_EMISSION
      LOGICAL, INTENT(IN)        :: DO_SURFACE_EMISSION

!  Deltam-2stream scaling flag

      LOGICAL, INTENT(IN)        :: DO_D2S_SCALING

!  BRDF surface flag

      LOGICAL, INTENT(IN)        :: DO_BRDF_SURFACE

!  Observational Geometry flag !@@ 2p1

      LOGICAL, INTENT(IN)        :: DO_USER_OBSGEOMS !@@ 2p1

!  @@ Rob Spurr, 23 January 2014, Version 2.3, SLEAVE option flags

      LOGICAL, INTENT(IN)        :: DO_SURFACE_LEAVING
      LOGICAL, INTENT(IN)        :: DO_SL_ISOTROPIC

!  BVP control --- New 6/25/14, Version 2.3 and higher
!   * PentaDiagonal Inverse flag (BVP solved from bottom to top). Only for BVPIndex = 1
!   * BVP Index : 0 = LAPACK, 1 = Penta # 1 (original), 2 = Penta # 2 (new, 2012 Kanal paper)
!   * BVP Scale Factor. Debug only. Set this to 1.0 on input

      LOGICAL      , INTENT(IN)  :: DO_PENTADIAG_INVERSE
      INTEGER      , INTENT(IN)  :: BVPINDEX
      REAL(kind=dp), INTENT(IN)  :: BVPSCALEFACTOR

!  4/12/24. CS THRESHOLD

      REAL(kind=dp), INTENT(IN)  :: CSTHRESH

!  Version 2.4, August 2014. Order of Taylor series (N) with Smallness number ( EPS)
!                            (Series including terms up to EPS^N)
      
      INTEGER      , intent(in)  :: TAYLOR_ORDER
      REAL(kind=dp), intent(in)  :: TAYLOR_SMALL

!  Thermal Cutoff (actually a layer optical thickness minimum)
!     Rob, introduced 14 May 2015, Version 2.4, following 2p3 implementation (2014)
!    Solutions are avoided for optically thin layers

      REAL(kind=dp), INTENT (IN) :: TCUTOFF

!  Numbers (basic), NTOTAL = 2 * NLAYERS

      INTEGER, INTENT(IN)        :: NLAYERS, NTOTAL

!  Stream value

      REAL(kind=dp), INTENT(IN)  :: STREAM_VALUE

!  Observational geometry input. [Same as LIDORT]. New 12/21/12 !@@ 2p1

      INTEGER, INTENT(IN)        :: N_USER_OBSGEOMS                    !@@ 2p1
      REAL(kind=dp), INTENT(IN)  :: USER_OBSGEOMS(MAX_USER_OBSGEOMS,3) !@@ 2p1

!  Viewing geometry. [Now Intent(inout), thanks to option for ObsGeom !@@ 2p1

      INTEGER, INTENT(INOUT)        :: N_USER_STREAMS
      REAL(kind=dp), INTENT(INOUT)  :: USER_ANGLES  ( MAX_USER_STREAMS )
      INTEGER, INTENT(INOUT)        :: N_USER_RELAZMS
      REAL(kind=dp), INTENT(INOUT)  :: USER_RELAZMS ( MAX_USER_RELAZMS )

!  Flux factor

      REAL(kind=dp), INTENT(IN)     :: FLUX_FACTOR

!  Solar geometry. [Now Intent(inout), thanks to option for ObsGeom !@@ 2p1

      INTEGER, INTENT(INOUT)        :: NBEAMS
      REAL(kind=dp), INTENT(INOUT)  :: BEAM_SZAS ( MAXBEAMS )

!  Height and earth radius (latter could be re-set internally)

      REAL(kind=dp), INTENT(INOUT) :: EARTH_RADIUS
      REAL(kind=dp), INTENT(IN)    :: HEIGHT_GRID ( 0:MAXLAYERS )

!  Geometry specification height
!      REAL(kind=dp), INTENT(IN)  :: GEOMETRY_SPECHEIGHT

!  Atmospheric optical properties

      REAL(kind=dp), INTENT(IN)  :: DELTAU_INPUT(MAXLAYERS)
      REAL(kind=dp), INTENT(IN)  :: OMEGA_INPUT (MAXLAYERS)
      REAL(kind=dp), INTENT(IN)  :: ASYMM_INPUT (MAXLAYERS)
      REAL(kind=dp), INTENT(IN)  :: D2S_SCALING (MAXLAYERS)

!  Atmospheric thermal sources

      REAL(kind=dp), INTENT(IN)  :: THERMAL_BB_INPUT ( 0:MAXLAYERS )

!  Surface properties

      REAL(kind=dp), INTENT(IN)  :: LAMBERTIAN_ALBEDO
      REAL(kind=dp), INTENT(IN)  :: SURFBB

!  BRDF fourier components
!  0 and 1 Fourier components of BRDF, following order (same all threads)
!    incident solar directions,  reflected quadrature stream
!    incident quadrature stream, reflected quadrature stream
!    incident solar directions,  reflected user streams    !  NOT REQUIRED
!    incident quadrature stream, reflected user streams

      REAL(kind=dp), INTENT(IN)  :: BRDF_F_0  ( 0:1, MAXBEAMS )
      REAL(kind=dp), INTENT(IN)  :: BRDF_F    ( 0:1 )
!      REAL(kind=dp), INTENT(IN)  :: UBRDF_F_0 ( 0:1, MAX_USER_STREAMS, MAXBEAMS )
      REAL(kind=dp), INTENT(IN)  :: UBRDF_F   ( 0:1, MAX_USER_STREAMS )

!  Surface thermal sources

      REAL(kind=dp), INTENT(IN)  :: EMISSIVITY

!  Version 2p3. 1/23/14. Introduce SLEAVE stuff
!  --------------------------------------------

!    Do not require any first-order inputs (exact or Fourier)

!  Isotropic Surface leaving term (if flag set)

      REAL(kind=dp), INTENT(IN) ::  SLTERM_ISOTROPIC ( MAXBEAMS )

!  Fourier components of Surface-leaving terms:
!    Every solar direction, SL-transmitted quadrature streams

      REAL(kind=dp), INTENT(IN) ::  SLTERM_F_0 ( 0:1, MAXBEAMS )

!  Exact Surface-Leaving term
!      REAL(kind=dp) ::  SLTERM_USERANGLES ( MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS )
!  Fourier components of Surface-leaving terms:
!    Every solar direction, SL-transmitted user streams. First order truncated
!      REAL(kind=dp) ::  USER_SLTERM_F_0 ( 0:1, MAX_USER_STREAMS, MAXBEAMS )

!  Output
!  ------

!  Radiance Results

      REAL(kind=dp), INTENT(INOUT) :: INTENSITY_TOA(MAX_GEOMETRIES)
      REAL(kind=dp), INTENT(INOUT) :: INTENSITY_BOA(MAX_GEOMETRIES)

!  Flux output
!     ! @@ Rob Spurr, 05 November 2013, Version 2.3 --> Flux Output

     REAL(kind=dp), INTENT(INOUT) :: FLUXES_TOA(MAXBEAMS,2)
     REAL(kind=dp), INTENT(INOUT) :: FLUXES_BOA(MAXBEAMS,2)

!  output solutions at ALL levels
!     ! @@ Rob Spurr, 17 July 2013, Version 2.2 --> Optional Output at ALL LEVELS

      REAL(kind=dp), INTENT(INOUT) :: RADLEVEL_UP (MAX_GEOMETRIES,0:MAXLAYERS)
      REAL(kind=dp), INTENT(INOUT) :: RADLEVEL_DN (MAX_GEOMETRIES,0:MAXLAYERS)

!  Fluxes at ALL levels (already initialized here)
!     Rob   Fix 5/24/17  for 2p4a --> Optional Output at ALL LEVELS

      REAL(kind=dp), INTENT(INOUT) :: FLUXLEVEL_DN (MAXBEAMS,0:MAXLAYERS,2)
      REAL(kind=dp), INTENT(INOUT) :: FLUXLEVEL_UP (MAXBEAMS,0:MAXLAYERS,2)

!  Numbers (geometry)
!   N_GEOMETRIES = NBEAMS * N_USER_STREAMS * N_USER_RELAZMS (Lattice value)
!   N_GEOMETRIES = N_USER_OBSGEOMS                          (OBsGeom value)

      INTEGER, INTENT(INOUT)       :: N_GEOMETRIES

!  Exception handling
!  ------------------

!    1. Check Messages and actions

      INTEGER      , INTENT(OUT) :: STATUS_INPUTCHECK
      INTEGER      , INTENT(OUT) :: C_NMESSAGES
      CHARACTER*100, INTENT(OUT) :: C_MESSAGES(0:MAXMESSAGES)
      CHARACTER*100, INTENT(OUT) :: C_ACTIONS (0:MAXMESSAGES)

!    2. Execution message and 2 Traces

      INTEGER      , INTENT(OUT) :: STATUS_EXECUTION
      CHARACTER*100, INTENT(OUT) :: E_MESSAGE, E_TRACE_1, E_TRACE_2

!  Local definitions
!  =================

!  Local Atmospheric Optical properties
!  ------------------------------------

!  After application of deltam scaling

      REAL(kind=dp) :: DELTAU_VERT(MAXLAYERS)
      REAL(kind=dp) :: OMEGA_TOTAL(MAXLAYERS)
      REAL(kind=dp) :: LOCAL_OMEGA_INPUT(MAXLAYERS)
      REAL(kind=dp) :: ASYMM_TOTAL(MAXLAYERS)

!  Chapman factors (from pseudo-spherical geometry)

      REAL(kind=dp) :: CHAPMAN_FACTORS ( MAXLAYERS, MAXLAYERS, MAXBEAMS )
      REAL(kind=dp) :: LOCAL_SZA       ( 0:MAXLAYERS, MAXBEAMS )

!  Miscsetup operations
!  ====================

!  Pseudo-spherical preparation
!  ----------------------------

!     Last layer to include Particular integral solution
!     Average-secant and initial tramsittance factors for solar beams.
!     Solar beam attenuation

      INTEGER       :: BEAM_CUTOFF ( MAXBEAMS )
      REAL(kind=dp) :: INITIAL_TRANS    ( MAXLAYERS, MAXBEAMS )
      REAL(kind=dp) :: AVERAGE_SECANT   ( MAXLAYERS, MAXBEAMS )
      REAL(kind=dp) :: LOCAL_CSZA       ( MAXLAYERS, MAXBEAMS )
      REAL(kind=dp) :: TRANS_SOLAR_BEAM( MAXBEAMS )

!  Derived optical thickness inputs

      REAL(kind=dp) :: DELTAU_SLANT ( MAXLAYERS, MAXLAYERS, MAXBEAMS )
      REAL(kind=dp) :: TAUSLANT     ( 0:MAXLAYERS, MAXBEAMS )

!  Reflectance flags

      LOGICAL       :: DO_DIRECTBEAM ( MAXBEAMS )

!  Transmittance Setups
!  --------------------

!  Transmittance factors for average secant stream

      REAL(kind=dp) :: T_DELT_MUBAR ( MAXLAYERS, MAXBEAMS )

!  Transmittance factors for user-defined stream angles
!    Computed in the initial setup stage for Fourier m = 0

      REAL(kind=dp) :: T_DELT_USERM ( MAXLAYERS, MAX_USER_STREAMS )
      REAL(kind=dp) :: ITRANS_USERM ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )

!  Forcing term multiplieres
!  -------------------------

!  coefficient functions for user-defined angles

      REAL(kind=dp) :: SIGMA_P(MAXLAYERS,MAX_USER_STREAMS,MAXBEAMS)
      REAL(kind=dp) :: SIGMA_M(MAXLAYERS,MAX_USER_STREAMS,MAXBEAMS)

!  L'Hopital's rule logical variables

      LOGICAL       :: EMULT_HOPRULE (MAXLAYERS,MAX_USER_STREAMS,MAXBEAMS)

!  Forcing term multipliers (saved for whole atmosphere)

      REAL(kind=dp) :: EMULT_UP (MAX_USER_STREAMS,MAXLAYERS,MAXBEAMS)
      REAL(kind=dp) :: EMULT_DN (MAX_USER_STREAMS,MAXLAYERS,MAXBEAMS)

!  Thermal help variables

      REAL(kind=dp) :: TCOM1       ( MAXLAYERS, 2 )
      REAL(kind=dp) :: THERMCOEFFS ( MAXLAYERS, 2 )

!  Fourier-component solutions
!  ===========================

      REAL(kind=dp) :: INTENSITY_F_UP (MAX_USER_STREAMS,MAXBEAMS)
      REAL(kind=dp) :: INTENSITY_F_DN (MAX_USER_STREAMS,MAXBEAMS)

!  Fourier-component solutions at ALL levels
!     ! @@ Rob Spurr, 17 July 2013, Version 2.2 --> Optional Output at ALL LEVELS

      REAL(kind=dp) :: RADLEVEL_F_UP (MAX_USER_STREAMS,MAXBEAMS,0:MAXLAYERS)
      REAL(kind=dp) :: RADLEVEL_F_DN (MAX_USER_STREAMS,MAXBEAMS,0:MAXLAYERS)

!  Single scatter solutions, commented out. THIS IS MS CODE !!!!
!      REAL(kind=dp) :: INTENSITY_SS_UP(N_GEOMETRIES)
!      REAL(kind=dp) :: INTENSITY_SS_DN(N_GEOMETRIES)

!  Other local variables
!  =====================

!  Local error handling

      CHARACTER(LEN=3) :: CF
      INTEGER          :: FOURIER, N_FOURIERS, STATUS_SUB
      INTEGER          :: N, UA, UM, IB, N_VIEWING, IBEAM, I, LUM, LUA, NCONS
      REAL(kind=dp)    :: AZM_ARGUMENT, DFC, DEG_TO_RAD, PI4
      REAL(kind=dp)    :: OMFAC, M1FAC, GDIFF, ALBEDO

!  Geometry offset arrays

      INTEGER          :: IBOFF ( MAXBEAMS )
      INTEGER          :: UMOFF ( MAXBEAMS, MAX_USER_STREAMS )

!  Post-processing flag (new for Version 2p3)

      LOGICAL          :: DO_POSTPROCESSING

!  Internal Full-quadrature control (true ==> Mu_bar = 1/sqrt(3) )

      LOGICAL          :: DO_FULLQUAD

!  Post-processing control mask

      INTEGER          :: N_PPSTREAMS, PPSTREAM_MASK ( MAX_USER_STREAMS, MAXBEAMS )

!  Local azimuth factors

      REAL(kind=dp)    :: AZMFAC (MAX_USER_STREAMS,MAXBEAMS,MAX_USER_RELAZMS)

!  Cosines and sines

      REAL(kind=dp)    :: X0  ( MAXBEAMS )
      REAL(kind=dp)    :: USER_STREAMS ( MAX_USER_STREAMS )
      REAL(kind=dp)    :: USER_SECANTS ( MAX_USER_STREAMS )
      REAL(kind=dp)    :: MUSTREAM, SINSTREAM

!  Input Conservative scattering flags

      LOGICAL          :: DO_CONSSCAT_INPUT(MAXLAYERS)

!mick - singularity buster output
      LOGICAL          :: SBUST(6)

!  Test variables

      !LOGICAL          :: DO_DEBUG_INPUT=.FALSE.
      !LOGICAL          :: DO_DEBUG_INPUT=.TRUE.

!  Initialize some variables
!  -------------------------

!  Input check

      STATUS_INPUTCHECK = 0
      C_NMESSAGES       = 0

!  Execution status and message/traces

      STATUS_EXECUTION  = 0
      E_MESSAGE = ' '
      E_TRACE_1 = ' '
      E_TRACE_2 = ' '

!  Local user indices; !@@ Only required for OBSGEOM option

      LUM = 1
      LUA = 1

!  Check input dimensions
!  ----------------------

      CALL TWOSTREAM_CHECK_INPUT_DIMS &
      ( DO_MVOUT_ONLY, DO_USER_OBSGEOMS, MAXMESSAGES, &
        MAXLAYERS, MAXTOTAL, MAXBEAMS, MAX_USER_STREAMS, MAX_USER_RELAZMS, MAX_USER_OBSGEOMS, &
        NLAYERS,   NTOTAL,   NBEAMS,   N_USER_STREAMS,   N_USER_RELAZMS,   N_USER_OBSGEOMS, &
        STATUS_SUB, C_NMESSAGES, C_MESSAGES, C_ACTIONS )

      IF ( STATUS_SUB .EQ. 1 ) THEN
        STATUS_INPUTCHECK = 1
        RETURN
      ENDIF

!  TWOSTREAM input debug
!  ---------------------

      IF (DO_DEBUG_INPUT) THEN
        CALL TWOSTREAM_DEBUG_INPUT_MASTER()
      END IF

!  Initialize output arrays
!  ------------------------

!mick fix 11/8/2012 - added main output
!  Main output

      INTENSITY_TOA(:) = zero
      INTENSITY_BOA(:) = zero

! @@ Rob Spurr, 17 July 2013, Version 2.2 --> Optional Output at ALL LEVELS

      RADLEVEL_UP (:,:) = zero
      RADLEVEL_DN (:,:) = zero

! @@ Rob Spurr, 05 November 2013, Version 2.3 --> BOA_TOA Flux outputs

      FLUXES_TOA(:,:) = zero
      FLUXES_BOA(:,:) = zero

!  Fluxes at ALL levels
!     Rob   Fix 5/24/17  for 2p4a --> Optional Output at ALL LEVELS

      FLUXLEVEL_UP(:,:,:) = zero
      FLUXLEVEL_DN(:,:,:) = zero

!  Constants
!  ---------

      DEG_TO_RAD = ACOS( - one ) / 180.0_dp
      PI4 = DEG_TO_RAD * 720.0_dp

!  Input checking
!  ==============

!  Check input Basic. This could be put outside the thread loop.
!    SS inputs are omitted in this version........
!    !@@ 2p1, Observational Geometry inputs are included (New 12/21/12)
!    !@@ 2p3 Includes check on Flux output flags, and setting of Post-Processing flag

!  5/8/24. Introduce flag DO_HEMISPHER_MEAN for thermal regime
!    -- Also stream_value (needed for this option). No deltam_scaling
!  5/9/24. Introduce flag DO_EDDINGTON_FLUX for solar flux regime
!    -- Also stream_value (needed for this option). No deltam_scaling
!  6/29/24. Insert flag DO_TRIDIAG, must be =  TRUE for DO_EDDINGTON_FLUX, DO_HEMISPHER_MEAN

      CALL TWOSTREAM_CHECK_INPUTS_BASIC  &
        ( MAXLAYERS, MAXMESSAGES, MAX_USER_OBSGEOMS,                     & ! Dimensions !@@
          MAXBEAMS, MAX_USER_STREAMS, MAX_USER_RELAZMS,                  & ! Dimensions
          DO_UPWELLING, DO_DNWELLING, DO_PLANE_PARALLEL, DO_D2S_SCALING, & ! Input
          DO_SOLAR_SOURCES, DO_THERMAL_EMISSION,                         & ! Input
          DO_HEMISPHER_MEAN, DO_EDDINGTON_FLUX, DO_TRIDIAG,              & ! Input (new options 5/9/24, 6/29/24)
          DO_MVOUT_ONLY, DO_ADDITIONAL_MVOUT, DO_POSTPROCESSING,         & ! Input !@@ New line, 2p3
          DO_USER_OBSGEOMS, N_USER_OBSGEOMS, USER_OBSGEOMS,              & ! Input !@@ New line
          NLAYERS, NBEAMS, N_USER_STREAMS, N_USER_RELAZMS,               & ! Input
          STREAM_VALUE, BEAM_SZAS, USER_ANGLES, USER_RELAZMS,            & ! Input
          EARTH_RADIUS, HEIGHT_GRID,                                     & ! Input
          STATUS_SUB, C_NMESSAGES, C_MESSAGES, C_ACTIONS )                 ! Output

      IF ( STATUS_SUB .EQ. 1 ) THEN
        STATUS_INPUTCHECK = 1
        RETURN
      ENDIF

!  Check input optical values (IOPs in the atmosphere)

     CALL TWOSTREAM_CHECK_INPUTS_OPTICAL &
       ( MAXLAYERS, MAXMESSAGES, NLAYERS, DO_TRIDIAG, DO_CLASSICAL_PI, & ! input
         DELTAU_INPUT, OMEGA_INPUT, ASYMM_INPUT,           & ! Input
         STATUS_SUB, C_NMESSAGES, C_MESSAGES, C_ACTIONS )    ! Output

      IF ( STATUS_SUB .EQ. 1 ) THEN
        STATUS_INPUTCHECK = 1
        RETURN
      ENDIF

!  Geometry offsets/masks
!  ======================

!  Save some offsets for indexing geometries

!   !@@ 2p1, This section revised for the Observational Geometry option
!   !@@ N_GEOMETRIES = NBEAMS * N_USER_STREAMS * N_USER_RELAZMS
!   !@@ 2p3, This section revised for the post-processing flag
!   !@@ 2p4, PPSTREAM masking for Lattice/Obsgeom choice

      N_VIEWING = 0 ; N_GEOMETRIES = 0
      IBOFF     = 0 ; UMOFF        = 0
      IF ( DO_USER_OBSGEOMS.and.DO_SOLAR_SOURCES ) THEN
         N_VIEWING    = N_USER_OBSGEOMS
         N_GEOMETRIES = N_USER_OBSGEOMS
      ELSE
         if ( DO_POSTPROCESSING ) THEN
            N_VIEWING    = N_USER_STREAMS * N_USER_RELAZMS
            N_GEOMETRIES = NBEAMS * N_VIEWING
            DO IBEAM = 1, NBEAMS
               IBOFF(IBEAM) = N_VIEWING * ( IBEAM - 1 )
               DO UM = 1, N_USER_STREAMS
                  UMOFF(IBEAM,UM) = IBOFF(IBEAM) +  N_USER_RELAZMS * (UM - 1)
               END DO
            END DO
         ENDIF
      ENDIF

!  Local post-processing control

      PPSTREAM_MASK = 0
      DO IB = 1, NBEAMS
         IF ( DO_USER_OBSGEOMS ) THEN
            N_PPSTREAMS = 1; PPSTREAM_MASK(1,IB) = IB
         else
            N_PPSTREAMS = N_USER_STREAMS
            do UM = 1, N_PPSTREAMS
               PPSTREAM_MASK(UM,IB) = UM
            enddo
         endif
      enddo

!  Geometry adjustment
!  -------------------

!  Not implemented. (only needed for Exact SS calculation)

!  Adjust surface condition
!      ADJUST_SURFACE = .FALSE.
!      IF ( DO_SSCORR_OUTGOING ) THEN
!        IF (HEIGHT_GRID(NLAYERS).GT.GEOMETRY_SPECHEIGHT ) THEN
!         ADJUST_SURFACE = .TRUE.
!        ENDIF
!      ENDIF
!  Perform adjustment.
!   Not implemented in streamlined version (only needed for SS stuff)
!      modified_eradius = earth_radius + GEOMETRY_SPECHEIGHT
!      CALL multi_outgoing_adjustgeom                                 &
!        ( N_USER_STREAMS, NBEAMS, N_USER_RELAZMS,                    &
!          N_USER_STREAMS,   NBEAMS,   N_USER_RELAZMS,                &
!          height_grid(nlayers), modified_eradius, adjust_surface,    &
!          user_angles,  beam_szas, user_relazms,                     &
!          user_angles_adjust, beam_szas_adjust, user_relazms_adjust, &
!          fail, mail )
!      if ( fail ) return

!  Chapman function calculation
!  ----------------------------

      DO IB = 1, NBEAMS
         CALL TWOSTREAM_BEAM_GEOMETRY_PREPARE &
            ( MAXLAYERS, MAXBEAMS,                      & ! Dimensions
              NLAYERS, DO_PLANE_PARALLEL, IB,           & ! Input
              BEAM_SZAS(IB), EARTH_RADIUS, HEIGHT_GRID, & ! Input
              CHAPMAN_FACTORS, LOCAL_SZA )                ! In/Out
      ENDDO

!  Get derived inputs
!  ==================

!  Quadrature. Add flag, 5/26/17

      DO_FULLQUAD = (ABS(STREAM_VALUE-0.5_dp).gt.1.0D-15)
      MUSTREAM  = STREAM_VALUE
      SINSTREAM = SQRT ( ONE - MUSTREAM * MUSTREAM )

!  Solar zenith angle cosine

      DO IB = 1, NBEAMS
         X0(IB) = COS ( BEAM_SZAS(IB) * DEG_TO_RAD )
      ENDDO

!  User stream cosines. 11/5/13 2p3 Post-processing control

      IF ( DO_POSTPROCESSING ) THEN
         DO I = 1, N_USER_STREAMS
            USER_STREAMS(I) = COS(DEG_TO_RAD * USER_ANGLES(I))
            USER_SECANTS(I) = ONE / USER_STREAMS(I)
         ENDDO
      ENDIF

!  Set Conservative scattering flags
!  4/12/24. Count the number of conservative-scatting layers
!  4/12/24. Introduce CSTHRESH = threshold for using Conservative scattering
!           Recommended value 1.0d-12 for solar, 1.0d-10 for thermal

!  6/29/24. Cannot do this for Eddington and Hemisphere mean cases
!           -- Otherwise Toggle for the quadrature cases

      DO_CONSSCAT_INPUT = .false. ; NCONS = 0
      if ( .not. DO_EDDINGTON_FLUX .and. .not. DO_HEMISPHER_MEAN ) then
        DO N = 1, NLAYERS 
!        if ( ABS(one-omega_input(n)) .lt. 1.0d-9 ) then
!        if ( ABS(one-omega_input(n)) .lt. 1.0d-15 ) then
          if ( ABS(one-omega_input(n)) .lt. CSTHRESH ) then
            DO_CONSSCAT_INPUT(n) = .true. ; NCONS = NCONS + 1
          endif
          local_omega_input(n) = omega_input(n)
        enddo
      ELSE
        DO N = 1, NLAYERS 
          if ( ABS(one-omega_input(n)) .lt. CSTHRESH ) then
            local_omega_input(n) = 1.0d0 - CSTHRESH
          else
            local_omega_input(n) = omega_input(n)
          endif
        ENDDO
     ENDIF

!  4/12/24. If all layers are conservative, and thermal emission only, then exit with warning
!  4/12/24. Works even with a surface emission term here

      if ( .not. DO_SOLAR_SOURCES .and. NCONS.eq.NLAYERS ) Then
        STATUS_EXECUTION = 2
        E_MESSAGE = 'Thermal Emission Only : All layers have conservative scattering!'
        E_TRACE_1 = 'Warning : 2S model still works but only with surface emission'
        E_TRACE_2 = 'Master program check, right after the DO_CONSSCAT_INPUT setup '
      endif

!  Set local atmospheric optical properties (Apply delta 2s scaling)
!  Just copy inputs, if not required
!   3/28/18. Deal with Conservative scattering option

      DO N = 1, NLAYERS
        IF ( DO_CONSSCAT_INPUT(n) ) THEN
          IF ( DO_D2S_SCALING ) THEN
            M1FAC = one - D2S_SCALING(N)
            GDIFF = ASYMM_INPUT(N) - D2S_SCALING(N)
            DELTAU_VERT(N) = M1FAC * DELTAU_INPUT(N)
            OMEGA_TOTAL(N) = ONE
            ASYMM_TOTAL(N) = GDIFF / M1FAC
          ELSE
            DELTAU_VERT(N) = DELTAU_INPUT(N)
            OMEGA_TOTAL(N) = ONE
            ASYMM_TOTAL(N) = ASYMM_INPUT(N)
          ENDIF
        ELSE
          IF ( DO_D2S_SCALING ) THEN
            OMFAC = one - OMEGA_INPUT(N) * D2S_SCALING(N)
            M1FAC = one - D2S_SCALING(N)
            GDIFF = ASYMM_INPUT(N) - D2S_SCALING(N)
            DELTAU_VERT(N) = OMFAC * DELTAU_INPUT(N)
            OMEGA_TOTAL(N) = M1FAC * local_omega_input(N) / OMFAC
            ASYMM_TOTAL(N) = GDIFF / M1FAC
          ELSE
            DELTAU_VERT(N) = DELTAU_INPUT(N)
            OMEGA_TOTAL(N) = local_omega_input(N)
            ASYMM_TOTAL(N) = ASYMM_INPUT(N)
          ENDIF
        ENDIF
      ENDDO

!mick fix 1/7/2012 - singularity busters added

!  Note: If running a case close to optical property numerical limits,
!        delta-m scaling may modify omega and/or g in such a way as to make
!        them unphysical or introduce instability; therefore, we recheck
!        omega and g AFTER delta-m scaling and slightly adjust them if necessary

      DO N = 1, NLAYERS
        SBUST = .false.
        go to 4444
        !Singularity buster for single scatter albedo
!        IF (OMEGA_TOTAL(N) > 0.999999999D0) THEN
!          OMEGA_TOTAL(N) = 0.999999999D0
!          SBUST(1) = .true.
!        ELSE IF (OMEGA_TOTAL(N) < 1.0D-9) THEN
!          OMEGA_TOTAL(N) = 1.0E-9_dp
!          SBUST(2) = .true.
!        END IF

        !Singularity buster for asymmetry parameter
        IF (ASYMM_TOTAL(N) > 0.999999999D0) THEN
          ASYMM_TOTAL(N) = 0.999999999D0
          SBUST(3) = .true.
        ELSE IF (ASYMM_TOTAL(N) < -0.999999999D0) THEN
          ASYMM_TOTAL(N) = -0.999999999D0
          SBUST(4) = .true.
        ELSE IF ((ASYMM_TOTAL(N) >= ZERO) .AND. &
                 (ASYMM_TOTAL(N) < 1.0D-9)) THEN
          ASYMM_TOTAL(N) = 1.0D-9
          SBUST(5) = .true.
        ELSE IF ((ASYMM_TOTAL(N) < ZERO) .AND. &
                 (ASYMM_TOTAL(N) > -1.0D-9)) THEN
          ASYMM_TOTAL(N) = -1.0D-9
          SBUST(6) = .true.
        END IF

        !WRITE(*,*)
        !WRITE(*,'(A,I2)') 'FOR LAYER: ',N
        !DO I=1,6
        !  WRITE(*,'(A,I1,A,L1)') '  SBUST(',I,') = ',SBUST(I)
        !ENDDO
4444  continue
      ENDDO

!  SETUP OPERATIONS (moved from Fourier, Version 2p3, 15 August 2014)
!  ==================================================================

!   MISCSETUPS (4 subroutines)  :
!       average-secant formulation,
!       transmittance setup
!       Thermal setup
!       Beam solution multipliers

!  Prepare quasi-spherical attenuation

      CALL TWOSTREAM_QSPREP &
        ( MAXLAYERS, MAXBEAMS, NLAYERS, NBEAMS, DO_PLANE_PARALLEL,     & ! Input
          DELTAU_VERT, CHAPMAN_FACTORS, X0, DO_DIRECTBEAM,             & ! In/Out
          INITIAL_TRANS, AVERAGE_SECANT, LOCAL_CSZA, BEAM_CUTOFF, & ! Output
          DELTAU_SLANT, TAUSLANT, TRANS_SOLAR_BEAM )                     ! Output

!  Transmittances and Transmittance factors. !@@ Add flag for Observation Geometry
!    !@@ Add Post-processing flag, 11/5/13

      CALL TWOSTREAM_PREPTRANS &
        ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS,                      & ! Dimensions
          DO_USER_OBSGEOMS, DO_POSTPROCESSING,                        & ! Input flags (2p1,2p3)
          NLAYERS, N_USER_STREAMS, NBEAMS, DELTAU_VERT, USER_SECANTS, & ! Input
          INITIAL_TRANS, AVERAGE_SECANT, BEAM_CUTOFF,            & ! Input
          T_DELT_MUBAR, T_DELT_USERM, ITRANS_USERM )                    ! Output

!  ** New, October 2011 **. Thermal setup
!  4/12/24. Add conservative scattering control

      IF ( DO_THERMAL_EMISSION ) THEN
         CALL TWOSTREAM_THERMALSETUP &
           ( MAXLAYERS, NLAYERS, DO_CONSSCAT_INPUT,      & ! Input
             OMEGA_TOTAL, DELTAU_VERT, THERMAL_BB_INPUT, & ! Input
             THERMCOEFFS, TCOM1 )                          ! input/Output
      ENDIF

!   EMULT_MASTER  : Beam source function multipliers.
!      !@@ Add alternative for Observational geometry, 2p1
!      !@@ Avoid altogether if no post-processing

!  Version 2.4 Overhaul-----
!     Rob  Fix 8/15/14  - Small numbers analysis using Taylor parameters
!     Rob  Fix 8/15/14  - Use of PPSTREAM and mask to deal with Obsgeom/Lattice choices
!     Rob  Fix 8/15/14  - Compact code in a single subroutine

      IF ( DO_SOLAR_SOURCES.and.DO_POSTPROCESSING ) THEN
         CALL TWOSTREAM_EMULTMASTER &
           ( MAXLAYERS, MAXBEAMS, MAX_USER_STREAMS,                  & ! Dimensions
             DO_UPWELLING, DO_DNWELLING, NLAYERS, NBEAMS,            & ! Input
             N_PPSTREAMS, PPSTREAM_MASK, TAYLOR_ORDER, TAYLOR_SMALL, & ! Input
             USER_SECANTS, DELTAU_VERT, T_DELT_MUBAR, T_DELT_USERM,  & ! Input
             BEAM_CUTOFF, ITRANS_USERM, AVERAGE_SECANT,         & ! Input
             SIGMA_M, SIGMA_P, EMULT_HOPRULE, EMULT_UP, EMULT_DN )     ! Output
      ENDIF

!      write(*,*)'EMULT',EMULT_UP(1,15,1), EMULT_DN(1,16,2)
!      write(*,*)'EMULT',EMULT_UP(2,23,2), EMULT_DN(2,1,1)

!  Initialise Fourier loop
!  =======================

!  Set Fourier number, Nominally 1 in absence of SS-only flag
!    Zero if no solar sources (Thermal-only run)
!  !@@ 2p3, Set NFOURIERS equal to zero for MVOUT_ONLY

      N_FOURIERS = 1
      IF (  DO_MVOUT_ONLY )         N_FOURIERS = 0
      IF ( .NOT. DO_SOLAR_SOURCES ) N_FOURIERS = 0

!mick fix 1/7/2012 - (test - make this permanent?)

      IF ( (NBEAMS == 1) .AND. (BEAM_SZAS(1) < 1.0D-8) ) N_FOURIERS = 0

!  Albedo

      ALBEDO = LAMBERTIAN_ALBEDO

!  Fourier loop
!  ============

      DO FOURIER = 0, N_FOURIERS

!  Azimuth cosine factors (Fourier = 1). !@@ 2p1, Notice OBSGEOM option
!  !@@ 2p3, not required for FLux-only output
!   5/8/24. Changed position of SOLAR_SOURCES flag

         AZMFAC = zero
         IF ( DO_POSTPROCESSING ) THEN
            IF ( FOURIER .GT. 0.and.DO_SOLAR_SOURCES ) THEN
               DFC = DBLE(FOURIER)
               IF ( DO_USER_OBSGEOMS ) THEN
                  DO IB = 1, NBEAMS
                     AZM_ARGUMENT = USER_RELAZMS(IB) * DFC
                     AZMFAC(LUM,IB,LUA) = COS(DEG_TO_RAD*AZM_ARGUMENT)
                  ENDDO
               ELSE
                  DO UA = 1, N_USER_RELAZMS
                     DO IB = 1, NBEAMS
                        DO UM = 1, N_USER_STREAMS
                           AZM_ARGUMENT = USER_RELAZMS(UA) * DFC
                           AZMFAC(UM,IB,UA) = COS(DEG_TO_RAD*AZM_ARGUMENT)
                        ENDDO
                     ENDDO
                  ENDDO
               ENDIF
            ENDIF
         ENDIF

!  Main call to Lidort Fourier module.
!  ----------------------------------

!  !@@ Add Observational Geometry dimension and control variables   !@@ 2p1
!  !@@ Call statement expanded to include ALL-LEVEL outputs         !@@ 2p2
!  !@@ Call statement expanded to include Flux outputs and control  !@@ 2p3  11/5/13
!  !@@ Call statement expanded to include Sleave inputs and control !@@ 2p3  01/23/14
!  !@@ Call statement expanded to include BVProblem control         !@@ 2p3  06/25/14
!  !@@ Call statement expanded to include Taylor-series control     !@@ 2p4  08/15/14
!  !@@ Call statement changed  to exclude Threading                 !@@ 2p4  08/15/14
!  !@@ Call statement changed  to include TCutoff                   !@@ 2p4  05/14/15

!  5/8/24. Introduce flag DO_HEMISPHER_MEAN for thermal regime
!   -- Thermal only, mubar = 0.5

!  5/9/24. Introduce flag DO_EDDINGTON_FLUX for solar regime, FLUXES-ONLY
!    -- Mubar must be 0.5, NNot allowed with thermal sources

         CALL TWOSTREAM_FOURIER_MASTER &
           ( MAXLAYERS, MAXTOTAL, MAXBEAMS, MAX_USER_STREAMS, DO_2S_LEVELOUT,   & ! Dimensions + flag
             DO_UPWELLING, DO_DNWELLING, DO_PLANE_PARALLEL, DO_BRDF_SURFACE,    & ! Input flags control
             DO_USER_OBSGEOMS, DO_CONSSCAT_INPUT, DO_SOLAR_SOURCES,             & ! Input flags control
             DO_THERMAL_EMISSION, DO_SURFACE_EMISSION,                          & ! Input flags
             DO_HEMISPHER_MEAN, DO_EDDINGTON_FLUX,                              & ! Input flags
             DO_MVOUT_ONLY, DO_ADDITIONAL_MVOUT, DO_POSTPROCESSING, DO_FULLQUAD,& ! Input flags !@@ New line, 2p3
             DO_SURFACE_LEAVING, DO_SL_ISOTROPIC, DO_PENTADIAG_INVERSE,         & ! Input flags !@@ 2p3 6/25/14
             BVPINDEX, BVPSCALEFACTOR, TAYLOR_ORDER, TAYLOR_SMALL, TCUTOFF,     & ! Input !@@ 2p3 6/25/14, 8/15/14
             NLAYERS, NTOTAL, NBEAMS, N_USER_STREAMS, FOURIER, PI4,             & ! Input integer control
             FLUX_FACTOR, STREAM_VALUE, X0, USER_STREAMS, USER_SECANTS,         & ! Input real control
             ALBEDO, BRDF_F_0, BRDF_F, UBRDF_F, DO_TRIDIAG, DO_CLASSICAL_PI,    & ! Input real surface
             SLTERM_ISOTROPIC, SLTERM_F_0, SURFBB, EMISSIVITY,                  & ! Input real sleave and thermal
             DELTAU_VERT, OMEGA_TOTAL, ASYMM_TOTAL, THERMCOEFFS,                & ! Input real optical
             BEAM_CUTOFF, INITIAL_TRANS, AVERAGE_SECANT, T_DELT_MUBAR,     & ! In/Out Miscsetups
             T_DELT_USERM, ITRANS_USERM, TRANS_SOLAR_BEAM, DO_DIRECTBEAM,       & ! In/Out Miscsetups
             SIGMA_P, SIGMA_M, EMULT_UP, EMULT_DN,                              & ! In/Out Miscsetups
             INTENSITY_F_UP, INTENSITY_F_DN, RADLEVEL_F_UP, RADLEVEL_F_DN,      & ! Output Radiances
             FLUXES_TOA, FLUXES_BOA, FLUXLEVEL_UP, FLUXLEVEL_DN,                & ! Output Fluxes
             STATUS_SUB, E_MESSAGE, E_TRACE_1 )                                   ! Exception handling

!  Exception handling

         IF ( STATUS_SUB .NE. 0 ) THEN
            STATUS_EXECUTION = 1
            WRITE(CF,'(I2)')FOURIER
            E_TRACE_2 = 'Error from 2S_FOURIER_MASTER, Fourier # ' //CF
            RETURN
         ENDIF

!  Fourier summation and Convergence examination
!  SS code not included in this version---------------
!  !@@ Alternative Call for Observationsl Geometry case      !@@ 2p1
!  !@@ Call statements expanded to include ALL-LEVEL outputs !@@ 2p2
!  !@@ Convergence skipped for MVOUT_ONLY option             !@@ 2p3 !mick fix 12/17/2013 - fixed logic

         IF ( .not. DO_MVOUT_ONLY ) then
            DO IBEAM = 1, NBEAMS
               IF ( DO_USER_OBSGEOMS.and.DO_SOLAR_SOURCES ) THEN
                  CALL TWOSTREAM_CONVERGE_OBSGEO &
                    ( MAXBEAMS, MAX_USER_STREAMS, MAX_USER_RELAZMS, & ! Dimensions
                      MAX_GEOMETRIES, MAXLAYERS,                    & ! Dimensions ! @@ 2p2
                      DO_UPWELLING, DO_DNWELLING, DO_2S_LEVELOUT,   & ! Inputs     ! @@ 2p2
                      NLAYERS, IBEAM, FOURIER, AZMFAC,              & ! Inputs     ! @@ 2p2
                      INTENSITY_F_UP,  INTENSITY_F_DN,              & ! Inputs
                      RADLEVEL_F_UP,   RADLEVEL_F_DN,               & ! Inputs     ! @@ 2p2
                      INTENSITY_TOA,   INTENSITY_BOA,               & ! In/Out
                      RADLEVEL_UP,     RADLEVEL_DN   )                ! In/Out     ! @@ 2p2
               ELSE
                  CALL TWOSTREAM_CONVERGE &
                    ( MAXBEAMS, MAX_USER_STREAMS, MAX_USER_RELAZMS,   & ! Dimensions
                      MAX_GEOMETRIES, MAXLAYERS,                      & ! Dimensions ! @@ 2p2
                      DO_UPWELLING, DO_DNWELLING, DO_2S_LEVELOUT,     & ! Inputs     ! @@ 2p2
                      NLAYERS, IBEAM, FOURIER,                        & ! Inputs     ! @@ 2p2
                      N_USER_STREAMS, N_USER_RELAZMS, AZMFAC, UMOFF,  & ! Inputs
                      INTENSITY_F_UP,  INTENSITY_F_DN,                & ! Inputs
                      RADLEVEL_F_UP,   RADLEVEL_F_DN,                 & ! Inputs     ! @@ 2p2
                      INTENSITY_TOA,   INTENSITY_BOA,                 & ! In/Out
                      RADLEVEL_UP,     RADLEVEL_DN   )                  ! In/Out     ! @@ 2p2
               ENDIF
            END DO
         ENDIF

!  End Fourier loop

      ENDDO

!  Finish

      RETURN

      CONTAINS

      SUBROUTINE TWOSTREAM_DEBUG_INPUT_MASTER()

!  May 2024. Add Hemisphere_mean and Eddington_Flux flags

      CALL TWOSTREAM_WRITE_STD_INPUT ( &
        MAXLAYERS, MAXTOTAL, MAXMESSAGES, MAXBEAMS, MAX_GEOMETRIES,     & 
        MAX_USER_RELAZMS, MAX_USER_STREAMS, MAX_USER_OBSGEOMS,          &
        DO_UPWELLING, DO_DNWELLING, DO_PLANE_PARALLEL, DO_2S_LEVELOUT,  &
        DO_MVOUT_ONLY, DO_ADDITIONAL_MVOUT, DO_CLASSICAL_PI,            &
        DO_SOLAR_SOURCES, DO_THERMAL_EMISSION, DO_SURFACE_EMISSION,     &
        DO_HEMISPHER_MEAN, DO_EDDINGTON_FLUX,                           & 
        DO_D2S_SCALING, DO_BRDF_SURFACE, DO_USER_OBSGEOMS,              &
        DO_SURFACE_LEAVING, DO_SL_ISOTROPIC, DO_PENTADIAG_INVERSE,      &
        BVPINDEX, BVPSCALEFACTOR, TAYLOR_ORDER, TAYLOR_SMALL, TCUTOFF,  &
        NLAYERS, NTOTAL, STREAM_VALUE, N_USER_OBSGEOMS, USER_OBSGEOMS,  &
        N_USER_STREAMS, USER_ANGLES, N_USER_RELAZMS, USER_RELAZMS,      &
        FLUX_FACTOR, NBEAMS, BEAM_SZAS, EARTH_RADIUS, HEIGHT_GRID,      &
        DELTAU_INPUT, OMEGA_INPUT, ASYMM_INPUT, D2S_SCALING,            &
        THERMAL_BB_INPUT, LAMBERTIAN_ALBEDO, SURFBB )

      IF (DO_BRDF_SURFACE) THEN
        CALL TWOSTREAM_WRITE_SUP_BRDF_INPUT (   &
          MAXBEAMS, MAX_USER_STREAMS, &
          NBEAMS, N_USER_STREAMS,     &
          BRDF_F_0, BRDF_F, UBRDF_F, EMISSIVITY )
      END IF

      IF (DO_SURFACE_LEAVING) THEN
        CALL TWOSTREAM_WRITE_SUP_SLEAVE_INPUT ( &
          MAXBEAMS,NBEAMS,&
          SLTERM_ISOTROPIC,SLTERM_F_0 )
      END IF

      END SUBROUTINE TWOSTREAM_DEBUG_INPUT_MASTER

END SUBROUTINE TWOSTREAM_MASTER

!

SUBROUTINE TWOSTREAM_FOURIER_MASTER &
        ( MAXLAYERS, MAXTOTAL, MAXBEAMS, MAX_USER_STREAMS, DO_2S_LEVELOUT,    & ! Dimensions + flag
          DO_UPWELLING, DO_DNWELLING, DO_PLANE_PARALLEL, DO_BRDF_SURFACE,     & ! Input flags control
          DO_USER_OBSGEOMS, DO_CONSSCAT_INPUT, DO_SOLAR_SOURCES,              & ! Input flags control
          DO_THERMAL_EMISSION, DO_SURFACE_EMISSION,                           & ! Input flags
          DO_HEMISPHER_MEAN, DO_EDDINGTON_FLUX,                               & ! Input flags
          DO_MVOUT_ONLY, DO_ADDITIONAL_MVOUT, DO_POSTPROCESSING, DO_FULLQUAD, & ! Input flags !@@ New line, 2p3
          DO_SURFACE_LEAVING, DO_SL_ISOTROPIC, DO_PENTADIAG_INVERSE,          & ! Input flags !@@ 2p3 6/25/14
          BVPINDEX, BVPSCALEFACTOR, TAYLOR_ORDER, TAYLOR_SMALL, TCUTOFF,      & ! Input !@@ 2p3 6/25/14, 8/15/14
          NLAYERS, NTOTAL, NBEAMS, N_USER_STREAMS, FOURIER, PI4,              & ! Input integer control
          FLUX_FACTOR, STREAM_VALUE, X0, USER_STREAMS, USER_SECANTS,          & ! Input real control
          ALBEDO, BRDF_F_0, BRDF_F, UBRDF_F, DO_TRIDIAG, DO_CLASSICAL_PI,     & ! Input real surface
          SLTERM_ISOTROPIC, SLTERM_F_0, SURFBB, EMISSIVITY,                   & ! Input real sleave and thermal
          DELTAU_VERT, OMEGA_TOTAL, ASYMM_TOTAL, THERMCOEFFS,                 & ! Input real optical
          BEAM_CUTOFF, INITIAL_TRANS, AVERAGE_SECANT, T_DELT_MUBAR,           & ! In/Out Miscsetups
          T_DELT_USERM, ITRANS_USERM, TRANS_SOLAR_BEAM, DO_DIRECTBEAM,        & ! In/Out Miscsetups
          SIGMA_P, SIGMA_M, EMULT_UP, EMULT_DN,                               & ! In/Out Miscsetups
          INTENSITY_F_UP, INTENSITY_F_DN, RADLEVEL_F_UP, RADLEVEL_F_DN,       & ! Output
          FLUXES_TOA, FLUXES_BOA, FLUXLEVEL_UP, FLUXLEVEL_DN,                 & ! Output Fluxes
          STATUS, MESSAGE, TRACE )                                              ! Exception handling

      implicit none

!  5/8/24. Introduce flag DO_HEMISPHER_MEAN for thermal regime
!   -- Thermal only, mubar = 0.5

!  5/9/24. Introduce flag DO_EDDINGTON_FLUX for solar regime, FLUXES-ONLY
!    -- Mubar must be 0.5, Not allowed with thermal sources

!  precision and parameters

      INTEGER      , PARAMETER :: dp   = KIND( 1.0D0 )
      REAL(kind=dp), parameter :: zero = 0.0_dp, one = 1.0_dp

!  input
!  -----

!  Dimensions :
!      MAXTOTAL  = 2 * MAXLAYERS
     
      INTEGER, INTENT(IN)        :: MAXLAYERS, MAXTOTAL
      INTEGER, INTENT(IN)        :: MAXBEAMS, MAX_USER_STREAMS

!  Flags
!     ! @@ Rob Spurr, 17 July 2013, Version 2.2, Levelout flag

      LOGICAL, INTENT(IN)  :: DO_UPWELLING, DO_DNWELLING
      LOGICAL, INTENT(IN)  :: DO_BRDF_SURFACE
      LOGICAL, INTENT(IN)  :: DO_2S_LEVELOUT

!  Plane parallel flag

      LOGICAL, INTENT(IN)  :: DO_PLANE_PARALLEL

!  Input Conservative scattering flags

      LOGICAL, INTENT(IN)  :: DO_CONSSCAT_INPUT(MAXLAYERS)

!  ** New **. October 2011, Sources control, including thermal

      LOGICAL, INTENT(IN)  :: DO_THERMAL_EMISSION
      LOGICAL, INTENT(IN)  :: DO_SURFACE_EMISSION
      LOGICAL, INTENT(IN)  :: DO_SOLAR_SOURCES

!  5/8/24. Introduce flag DO_HEMISPHER_MEAN for thermal regime

      Logical, intent(in) :: DO_HEMISPHER_MEAN

!  5/9/24. Introduce flag DO_EDDINGTON_FLUX for solar flux regime

      Logical, intent(in) :: DO_EDDINGTON_FLUX

!   !@@ Observational Geometry flag !@@ 2p1

      LOGICAL, INTENT(IN)  :: DO_USER_OBSGEOMS !@@ 2p1

!  !@@ Version 2p3, 11/5/13. Flux output flags, processing flag

      LOGICAL, INTENT(IN)  :: DO_MVOUT_ONLY
      LOGICAL, INTENT(IN)  :: DO_ADDITIONAL_MVOUT
      LOGICAL, INTENT(IN)  :: DO_POSTPROCESSING

!   Full quadrature control (true ==> Mu_bar = 1/sqrt(3) )

      LOGICAL, INTENT(IN)  :: DO_FULLQUAD

!  !@@ Version 2p3, 1/23/14. Surface leaving control

      LOGICAL, INTENT(IN)  :: DO_SURFACE_LEAVING
      LOGICAL, INTENT(IN)  :: DO_SL_ISOTROPIC

!  BVP control --- New 6/25/14, Version 2.3 and higher
!   * PentaDiagonal Inverse flag (BVP solved from bottom to top). Only for BVPIndex = 1
!   * BVP Index : 0 = LAPACK, 1 = Penta # 1 (original), 2 = Penta # 2 (new, 2012 Kanal paper)
!   * BVP Scale Factor. Debug only. Set this to 1.0 on input

      LOGICAL      , INTENT(IN)  :: DO_PENTADIAG_INVERSE
      INTEGER      , INTENT(IN)  :: BVPINDEX
      REAL(kind=dp), INTENT(IN)  :: BVPSCALEFACTOR

!  Order of Taylor series (including terms up to EPS^n)
      
      INTEGER, intent(in)        :: TAYLOR_ORDER
      REAL(kind=dp), intent(in)  :: TAYLOR_SMALL

!  Thermal Cutoff (actually a layer optical thickness minimum)
!     Rob, introduced 14 May 2015, following 2p3 implementation (2014)
!    Solutions are avoided for optically thin layers

      REAL(kind=dp), INTENT (IN) :: TCUTOFF

!  Numbers

      INTEGER, INTENT(IN)  :: NLAYERS, NTOTAL
      INTEGER, INTENT(IN)  :: NBEAMS, N_USER_STREAMS

!  Input Fourier component number

      INTEGER, INTENT(IN)        :: FOURIER

!  4pi

      REAL(kind=dp), INTENT(IN)  :: PI4

!  Flux factor

      REAL(kind=dp), INTENT(IN)  :: FLUX_FACTOR

!  Stream value

      REAL(kind=dp), INTENT(IN)  :: STREAM_VALUE

!  Geometry

      REAL(kind=dp), INTENT(IN)  :: X0           ( MAXBEAMS )
      REAL(kind=dp), INTENT(IN)  :: USER_STREAMS ( MAX_USER_STREAMS )
      REAL(kind=dp), INTENT(IN)  :: USER_SECANTS ( MAX_USER_STREAMS )

!  Surface variables
!  ------------------

!  Lambertian Albedo

     REAL(kind=dp), INTENT(IN)  :: ALBEDO

!  BRDF Fourier components (NOT threaded)
!  0 and 1 Fourier components of BRDF, following order (same all threads)
!    incident solar directions,  reflected quadrature stream
!    incident quadrature stream, reflected quadrature stream
!    incident solar directions,  reflected user streams -- NOT REQUIRED
!    incident quadrature stream, reflected user streams

      REAL(kind=dp), INTENT(IN)  :: BRDF_F_0  ( 0:1, MAXBEAMS )
      REAL(kind=dp), INTENT(IN)  :: BRDF_F    ( 0:1 )
!      REAL(kind=dp), INTENT(IN)  :: UBRDF_F_0 ( 0:1, MAX_USER_STREAMS, MAXBEAMS )
      REAL(kind=dp), INTENT(IN)  :: UBRDF_F   ( 0:1, MAX_USER_STREAMS )

!  New Tridiagonal flag

      LOGICAL, INTENT(IN)        :: DO_TRIDIAG

!  Classical PI flag
!  Rob Fix 5/24/17 for 2p4a, reintroduce classical PI methods

      LOGICAL, INTENT(IN)        :: DO_CLASSICAL_PI

!  ** New **. October 2011. Thermal variables
!  ------------------------------------------

      REAL(kind=dp), INTENT(IN)  :: SURFBB
      REAL(kind=dp), INTENT(IN)  :: EMISSIVITY

!  Version 2p3. 1/23/14. Introduce SLEAVE stuff
!  --------------------------------------------

!    Do not require any first-order inputs (exact or Fourier)

!  Isotropic Surface leaving term (if flag set)

      REAL(kind=dp), INTENT(IN) ::  SLTERM_ISOTROPIC ( MAXBEAMS )

!  Fourier components of Surface-leaving terms:
!    Every solar direction, SL-transmitted quadrature streams

      REAL(kind=dp), INTENT(IN) ::  SLTERM_F_0 ( 0:1, MAXBEAMS )

!  Optical properties
!  ------------------

      REAL(kind=dp), INTENT(IN)  :: DELTAU_VERT(MAXLAYERS)
      REAL(kind=dp), INTENT(IN)  :: OMEGA_TOTAL(MAXLAYERS)
      REAL(kind=dp), INTENT(IN)  :: ASYMM_TOTAL(MAXLAYERS)

!  Output
!  ------

!  User-defined solutions

      REAL(kind=dp), INTENT(OUT) :: INTENSITY_F_UP (MAX_USER_STREAMS,MAXBEAMS)
      REAL(kind=dp), INTENT(OUT) :: INTENSITY_F_DN (MAX_USER_STREAMS,MAXBEAMS)

!  Flux output (already initialized here)
!     ! @@ Rob Spurr, 05 November 2013, Version 2.3 --> Flux Output

     REAL(kind=dp), INTENT(INOUT) :: FLUXES_TOA(MAXBEAMS,2)
     REAL(kind=dp), INTENT(INOUT) :: FLUXES_BOA(MAXBEAMS,2)

!  Fourier-component solutions at ALL levels
!     ! @@ Rob Spurr, 17 July 2013, Version 2.2 --> Optional Output at ALL LEVELS

      REAL(kind=dp), INTENT(OUT) :: RADLEVEL_F_UP (MAX_USER_STREAMS,MAXBEAMS,0:MAXLAYERS)
      REAL(kind=dp), INTENT(OUT) :: RADLEVEL_F_DN (MAX_USER_STREAMS,MAXBEAMS,0:MAXLAYERS)

!  Fluxes at ALL levels (already initialized here)
!     Rob   Fix 5/24/17  for 2p4a --> Optional Output at ALL LEVELS

      REAL(kind=dp), INTENT(INOUT) :: FLUXLEVEL_DN (MAXBEAMS,0:MAXLAYERS,2)
      REAL(kind=dp), INTENT(INOUT) :: FLUXLEVEL_UP (MAXBEAMS,0:MAXLAYERS,2)

!  Single scatter solutions, commented out in this streamlined version
!      REAL(kind=dp) INTENSITY_SS_UP (N_GEOMETRIES)
!      REAL(kind=dp) INTENSITY_SS_DN (N_GEOMETRIES)

!  Exception handling

      INTEGER      , INTENT(OUT)  :: STATUS
      CHARACTER*(*), INTENT(OUT)  :: MESSAGE, TRACE

!  Miscsetups Arrays required as In/OUT
!  ====================================

!  Thermal help variables
!  ----------------------

      REAL(kind=dp), INTENT(INOUT) :: THERMCOEFFS ( MAXLAYERS, 2 )

!  Solar beam pseudo-spherical setup
!  ---------------------------------

!     Last layer to include Particular integral solution
!     Average-secant and initial/layer tramsittance factors for solar beams.

      INTEGER      , INTENT(INOUT) :: BEAM_CUTOFF      ( MAXBEAMS )
      REAL(kind=dp), INTENT(INOUT) :: INITIAL_TRANS    ( MAXLAYERS, MAXBEAMS )
      REAL(kind=dp), INTENT(INOUT) :: AVERAGE_SECANT   ( MAXLAYERS, MAXBEAMS )
      REAL(kind=dp), INTENT(INOUT) :: T_DELT_MUBAR     ( MAXLAYERS, MAXBEAMS )

!  Solar beam attenuation, reflectance flag
!  ----------------------------------------

      REAL(kind=dp), INTENT(INOUT) :: TRANS_SOLAR_BEAM ( MAXBEAMS )
      LOGICAL      , INTENT(INOUT) :: DO_DIRECTBEAM ( MAXBEAMS )

!  Transmittance for user-defined stream angles
!  --------------------------------------------

      REAL(kind=dp), INTENT(INOUT) :: ITRANS_USERM ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )
      REAL(kind=dp), INTENT(INOUT) :: T_DELT_USERM ( MAXLAYERS, MAX_USER_STREAMS )

!  Multiplier arrays
!  -----------------

!  coefficient functions for user-defined angles

      REAL(kind=dp), INTENT(INOUT) :: SIGMA_P(MAXLAYERS,MAX_USER_STREAMS,MAXBEAMS)
      REAL(kind=dp), INTENT(INOUT) :: SIGMA_M(MAXLAYERS,MAX_USER_STREAMS,MAXBEAMS)

!  Forcing term multipliers (saved for whole atmosphere)

      REAL(kind=dp), INTENT(INOUT) :: EMULT_UP (MAX_USER_STREAMS,MAXLAYERS,MAXBEAMS)
      REAL(kind=dp), INTENT(INOUT) :: EMULT_DN (MAX_USER_STREAMS,MAXLAYERS,MAXBEAMS)

!  Local Arrays for Use in Subroutines
!  ===================================

!  Geometry arrays
!  ---------------

!  These just save some Polynomial expansions

      REAL(kind=dp) :: ULP  ( MAX_USER_STREAMS )
      REAL(kind=dp) :: POX  ( MAXBEAMS )
      REAL(kind=dp) :: PX0X ( MAXBEAMS )
      REAL(kind=dp) :: PX11, PXSQ

!  Solar beam Attenuation
!  ----------------------

!  Atmospheric attenuation

      REAL(kind=dp) :: ATMOS_ATTN ( MAXBEAMS )

!  Direct beam solutions. No USER-term required, MS-mode only

      REAL(kind=dp) :: DIRECT_BEAM ( MAXBEAMS )

!  Multiplier arrays (Homogeneous solutions)
!  -----------------

!  coefficient functions for user-defined angles

      REAL(kind=dp) :: ZETA_M(MAX_USER_STREAMS,MAXLAYERS)
      REAL(kind=dp) :: ZETA_P(MAX_USER_STREAMS,MAXLAYERS)

!  Integrated homogeneous solution multipliers, whole layer

      REAL(kind=dp) :: HMULT_1(MAX_USER_STREAMS,MAXLAYERS)
      REAL(kind=dp) :: HMULT_2(MAX_USER_STREAMS,MAXLAYERS)

!  Solutions to the homogeneous RT equations
!  -----------------------------------------

!  local matrices for eigenvalue computation

      REAL(kind=dp) :: SAB(MAXLAYERS), DAB(MAXLAYERS)

!  Gamma, Alpha and Epsilon

      REAL(kind=dp)  :: GAMMA(MAXLAYERS)
      REAL(kind=dp)  :: ALPHA(MAXLAYERS)
      REAL(kind=dp)  :: EPSON(MAXLAYERS)

!  Eigensolutions

      REAL(kind=dp) :: EIGENVALUE(MAXLAYERS)
      REAL(kind=dp) :: EIGENTRANS(MAXLAYERS)

!  Eigenvector solutions

      REAL(kind=dp) :: XPOS(2,MAXLAYERS)
      REAL(kind=dp) :: XNEG(2,MAXLAYERS)

!  Green;s function normalization factors
!    Introduced for [V2p3, Mark 10]

      REAL(kind=dp) :: NORM_SAVED(MAXLAYERS)

!  Specialist solutions for DO_HEMISPHER_MEAN or DO_EDDINGTON_FLUX

      REAL(kind=dp) :: GAMMA1(MAXLAYERS)
      REAL(kind=dp) :: GAMMA2(MAXLAYERS)
      REAL(kind=dp) :: BIG_GAMMA(MAXLAYERS)

!  Saved help variables

      REAL(kind=dp) :: U_HELP_P(0:1)
      REAL(kind=dp) :: U_HELP_M(0:1)

!  Eigenvectors defined at user-defined stream angles
!     EP for the positive KEIGEN values, EM for -ve KEIGEN

      REAL(kind=dp) :: U_XPOS(MAX_USER_STREAMS,MAXLAYERS)
      REAL(kind=dp) :: U_XNEG(MAX_USER_STREAMS,MAXLAYERS)

!  Downwelling BOA solution, before reflectance

      REAL(kind=dp) :: H_HOMP, H_HOMM

!  Boundary Value Problem
!  Original and Elimination matrices (Pentadiagonal, 2x2)

      REAL(kind=dp) :: SELM   (2,2)
      REAL(kind=dp) :: ELM    (MAXTOTAL,4)
      REAL(kind=dp) :: MAT    (MAXTOTAL,5)

!  particular integrals and BVP solution
!  -------------------------------------

!  Solutions at layer boundaries

      REAL(kind=dp) :: WUPPER(2,MAXLAYERS)
      REAL(kind=dp) :: WLOWER(2,MAXLAYERS)

!  Downwelling BOA solution, before reflectance

      REAL(kind=dp) :: H_PARTIC

!  Single-scatter Particular beam solutions at user-defined angles
!    ****** NOT REQUIRED for MS-mode only
!      REAL(kind=dp) :: U_WPOS1(MAX_USER_STREAMS,MAXLAYERS)
!      REAL(kind=dp) :: U_WNEG1(MAX_USER_STREAMS,MAXLAYERS)

!  Solution constants of integration, and related quantities

      REAL(kind=dp) :: LCON(MAXLAYERS)
      REAL(kind=dp) :: MCON(MAXLAYERS)
      REAL(kind=dp) :: LCON_XVEC(2,MAXLAYERS)
      REAL(kind=dp) :: MCON_XVEC(2,MAXLAYERS)

!  Beam Solutions (Classical)
!  --------------------------

!  Rob Fix 5/24/17 for 2p4a, reintroduce classical PI methods

!  Auxiliary vectors

      REAL(kind=dp)  :: QDIFVEC(MAXLAYERS)
      REAL(kind=dp)  :: QSUMVEC(MAXLAYERS)
      REAL(kind=dp)  :: QVEC   (MAXLAYERS)

!  Beam solution

      REAL(kind=dp)  :: WVEC(2,MAXLAYERS)

!  Diffuse-term Particular beam solutions at user-defined angles

      REAL(kind=dp) :: U_WPOS2(N_USER_STREAMS,NLAYERS)
      REAL(kind=dp) :: U_WNEG2(N_USER_STREAMS,NLAYERS)

!  Saved help variables

      REAL(kind=dp) :: W_HELP(0:1)

!  Beam Solutions (Greens function)
!  --------------------------------

!  Saved quantities for the Green function solution

      REAL(kind=dp) :: ATERM_SAVE(MAXLAYERS)
      REAL(kind=dp) :: BTERM_SAVE(MAXLAYERS)
      REAL(kind=dp) :: DMI, DPI

!  Layer C and D functions

      REAL(kind=dp) :: CFUNC(MAXLAYERS)
      REAL(kind=dp) :: DFUNC(MAXLAYERS)

!  Green function Multipliers for solution

      REAL(kind=dp) :: GFUNC_UP(MAXLAYERS)
      REAL(kind=dp) :: GFUNC_DN(MAXLAYERS)

!  Holding arrays for Multiplier coefficients

      REAL(kind=dp) :: GAMMA_M(MAXLAYERS)
      REAL(kind=dp) :: GAMMA_P(MAXLAYERS)

!  Thermal solutions
!  -----------------

!  Saved quantities for the Green function solution

      REAL(kind=dp) :: TTERM_SAVE (MAXLAYERS)
      REAL(kind=dp) :: T_C_MINUS (MAXLAYERS,0:2)
      REAL(kind=dp) :: T_C_PLUS  (MAXLAYERS,0:2)

!  Thermal solution at layer boundaries

      REAL(kind=dp) :: T_WUPPER ( 2, MAXLAYERS )
      REAL(kind=dp) :: T_WLOWER ( 2, MAXLAYERS )

!  post-processed thermal solutions (Classical Pi only). Rob Fix 5/24/17

      REAL(kind=dp) :: U_TPOS2 ( MAX_USER_STREAMS, MAXLAYERS, 2 )
      REAL(kind=dp) :: U_TNEG2 ( MAX_USER_STREAMS, MAXLAYERS, 2 )

!  Complete layer term solutions

      REAL(kind=dp) :: LAYER_TSUP_UP(MAX_USER_STREAMS,MAXLAYERS)
      REAL(kind=dp) :: LAYER_TSUP_DN(MAX_USER_STREAMS,MAXLAYERS)

!  Post-processing variables
!  -------------------------

!  Greens function multipliers
!  Rob Fix 8/15/14. Source function integrated Green function multipliers (whole layer)

      REAL(kind=dp) :: PMULT_DU(MAX_USER_STREAMS,MAXLAYERS)
      REAL(kind=dp) :: PMULT_DD(MAX_USER_STREAMS,MAXLAYERS)
      REAL(kind=dp) :: PMULT_UU(MAX_USER_STREAMS,MAXLAYERS)
      REAL(kind=dp) :: PMULT_UD(MAX_USER_STREAMS,MAXLAYERS)

!  Reflectance integrand  a(j).x(j).I(-j)

      REAL(kind=dp) :: IDOWNSURF

!  Cumulative source terms

      REAL(kind=dp) :: CUMSOURCE_UP(MAX_USER_STREAMS,0:MAXLAYERS)
      REAL(kind=dp) :: CUMSOURCE_DN(MAX_USER_STREAMS,0:MAXLAYERS)

!  Local help variables
!  --------------------

      INTEGER :: N, IBEAM, i

!  local inclusion flags. ** New October 2011 **, thermal flags
! !@@ 2p3 11/5/13. Control for the Flux calculation  

      LOGICAL :: DO_INCLUDE_SURFACE
      LOGICAL :: DO_INCLUDE_SURFEMISS
      LOGICAL :: DO_INCLUDE_THERMEMISS
      LOGICAL :: DO_INCLUDE_DIRECTBEAM
      LOGICAL :: DO_INCLUDE_MVOUT

!  Flux multiplier and Fourier component numbers

      REAL(kind=dp) :: FLUXMULT
      REAL(kind=dp) :: DELTA_FACTOR
      REAL(kind=dp) :: SURFACE_FACTOR

!  Error tracing

      INTEGER       :: STATUS_SUB

!  4/22/17. Tridiagonal tester
!  ---------------------------

!  control flag

      !LOGICAL          :: DO_TRIDIAG = .FALSE.
      !LOGICAL          :: DO_TRIDIAG = .TRUE.

!  Tridiagonal Matrix entries for solving BCs

      REAL(kind=dp) :: AMAT_TD(MAXTOTAL)
      REAL(kind=dp) :: BMAT_TD(MAXTOTAL)
      REAL(kind=dp) :: DMAT_TD(MAXTOTAL)

!  Tridiagonal help matrices

      REAL(kind=dp) :: EMAT_TD (MAXLAYERS,4)

!  Conservative scattering solutions and flags. Postprocessed sources, Introduced 3/28/18

      REAL(kind=dp) :: TGMFUNC(MAXLAYERS),TGM(MAXLAYERS)
      LOGICAL       :: CONSSCAT(MAXLAYERS)
      REAL(kind=dp) :: CSPPHOMUP(2,MAX_USER_STREAMS,MAXLAYERS)
      REAL(kind=dp) :: CSPPHOMDN(2,MAX_USER_STREAMS,MAXLAYERS)

!  ##############
!  initialization
!  ##############

!  Exception handling initialize

      STATUS = 0
      MESSAGE = ' '
      TRACE   = ' '

!  Initialize new output. NOT EFFICIENT - MICK, any suggestions ???
!     ! @@ Rob Spurr, 17 July 2013, Version 2.2 --> Optional output at ALL LEVELS

      RADLEVEL_F_UP = zero
      RADLEVEL_F_DN = zero

!  Set local flags
!  ---------------

!  inclusion of thermal surface emission term, only for Fourier = 0

      DO_INCLUDE_SURFEMISS = .FALSE.
      IF ( DO_SURFACE_EMISSION ) THEN
        IF ( FOURIER .EQ. 0 ) THEN
          DO_INCLUDE_SURFEMISS = .TRUE.
        ENDIF
      ENDIF

!  inclusion of thermal emission term, only for Fourier = 0

      DO_INCLUDE_THERMEMISS = .FALSE.
      IF ( DO_THERMAL_EMISSION ) THEN
        IF ( FOURIER .EQ. 0 ) THEN
          DO_INCLUDE_THERMEMISS = .TRUE.
        ENDIF
      ENDIF

!  Surface flag (for inclusion of some kind of reflecting boundary)

      DO_INCLUDE_SURFACE = .FALSE.
      IF ( DO_BRDF_SURFACE ) THEN
        DO_INCLUDE_SURFACE = .TRUE.
      ELSE
!mick fix 1/30/2015 - refined control logic
        IF ( FOURIER .EQ. 0 .AND. ALBEDO .NE. ZERO) DO_INCLUDE_SURFACE = .TRUE.
        !IF ( FOURIER .EQ. 0 ) DO_INCLUDE_SURFACE = .TRUE.
      ENDIF

!  Direct beam flag (only if above surface flag has been set)

!mick fix 1/30/2015 - refined control logic
      IF ( DO_SOLAR_SOURCES .and. DO_INCLUDE_SURFACE ) THEN
      !IF ( DO_INCLUDE_SURFACE ) THEN
        DO IBEAM = 1, NBEAMS
          DO_DIRECTBEAM(IBEAM) = .TRUE.
        ENDDO
      ELSE
        DO IBEAM = 1, NBEAMS
          DO_DIRECTBEAM(IBEAM) = .FALSE.
        ENDDO
      ENDIF

!  Inclusion of mean value calculation
! !@@ 2p3 11/5/13. Control for the Flux calculation  

      DO_INCLUDE_MVOUT = .FALSE.
      IF ( DO_ADDITIONAL_MVOUT .OR. DO_MVOUT_ONLY ) THEN
        IF ( FOURIER .EQ. 0 ) THEN
          DO_INCLUDE_MVOUT = .TRUE.
        ENDIF
      ENDIF

!  surface reflectance factors

      IF ( FOURIER .EQ. 0 ) THEN
        SURFACE_FACTOR = 2.0_dp
        DELTA_FACTOR   = one
      ELSE
        SURFACE_FACTOR = one
        DELTA_FACTOR   = 2.0_dp
      ENDIF

!  Flux multiplier
!   = 1 / 4.pi with beam sources, 1.0 for thermal

      FLUXMULT   = DELTA_FACTOR

!  Conservative scattering. Only for M = 0
!  4/12/14. More control here on CS terms necessary in layers with conservative scattering.

      CONSSCAT = .false. ; TGM = zero ; TGMFUNC = zero
      IF ( FOURIER .EQ. 0 ) THEN
        CONSSCAT(1:NLAYERS) = DO_CONSSCAT_INPUT(1:NLAYERS)
        DO N = 1, NLAYERS 
          IF ( CONSSCAT(N) ) then
            TGMFUNC(N) = ( one - 3.d0*STREAM_VALUE**2*ASYMM_TOTAL(N) ) / 2.d0 / STREAM_VALUE
            TGM(N) = deltau_VERT(n) * TGMFUNC(N)
          ENDIF
        ENDDO
      ENDIF

!  Old code
!      CONSSCAT = .false. ; TGM = zero ; TGMFUNC = zero
!      if ( fourier.eq.0 ) then
!         DO N = 1, NLAYERS 
!!            if ( ABS(one-omega_total(n)) .lt. 1.0d-9 ) then
!            if ( ABS(one-omega_total(n)) .lt. 1.0d-13 ) then
!!            if ( ABS(one-omega_total(n)) .lt. 1.0d-15 ) then
!               CONSSCAT(n) = .true.
!               TGMFUNC(N) = ( one - 3.d0*STREAM_VALUE**2*ASYMM_TOTAL(N) ) / 2.d0 / STREAM_VALUE
!               TGM(N) = deltau_VERT(n) * TGMFUNC(N)
!            endif
!         enddo
!      endif

!  Reflected Direct beam attenuation.
!  ! @@2p3, 1/23/14 add SLEAVE inputs
!  6/29/24. Include Eddington Flux flag in the argument list

      CALL TWOSTREAM_DIRECTBEAM & 
        ( MAXBEAMS, DO_EDDINGTON_FLUX,          & ! Dimensions
          DO_INCLUDE_SURFACE, DO_BRDF_SURFACE,  & ! Input
          DO_SURFACE_LEAVING, DO_SL_ISOTROPIC,  & ! input @@ 2p3
          NBEAMS, FOURIER, FLUX_FACTOR, X0, STREAM_VALUE, & ! Input
          DELTA_FACTOR, ALBEDO, BRDF_F_0,       & ! Input
          SLTERM_ISOTROPIC, SLTERM_F_0,         & ! input @@ 2p3
          TRANS_SOLAR_BEAM, DO_DIRECTBEAM,      & ! Input
          ATMOS_ATTN, DIRECT_BEAM )               ! Output

!  Auxiliary Geometry
!  ! @@2p3, 11/5/13 add Post-processing flag
!   5/9/24. Will not be done for the FLux-only alternative methods (HemisphericMean, Eddington)

      if ( .not. ( DO_HEMISPHER_MEAN .or. DO_EDDINGTON_FLUX) ) THEN
        CALL TWOSTREAM_AUXGEOM &
          ( MAX_USER_STREAMS, MAXBEAMS, DO_POSTPROCESSING, & ! Dimensions, Flag (2p3
            N_USER_STREAMS, NBEAMS, FOURIER, & ! Input
            X0, USER_STREAMS, STREAM_VALUE,  & ! Input
            PX11, PXSQ, POX, PX0X, ULP )       ! Output
      endif

!  #########################################
!   RTE HOMOGENEOUS SOLUTIONS and BVP SETUP
!  #########################################

!  Start layer loop

      DO N = 1, NLAYERS

!  Get Discrete ordinate solutions for this layer
!    Version 2p4. Green's function output = NORM_SAVED
!  4/12/24. Zero the solutions if there is conservative scattering in the layer.
!  5/9/24. Introduce alternative HOM_SOLUTION for 2 other flux options
!  7/16/24. BIG_GAMMA = GAMMA

         if ( DO_HEMISPHER_MEAN .or. DO_EDDINGTON_FLUX) THEN
           CALL OTHER2S_HOM_SOLUTION &
          ( MAXLAYERS, DO_HEMISPHER_MEAN, DO_EDDINGTON_FLUX, & ! Inputs
            N, OMEGA_TOTAL, ASYMM_TOTAL, DELTAU_VERT,          & ! Inputs
            GAMMA1, GAMMA2, BIG_GAMMA, EIGENVALUE, EIGENTRANS, XPOS, XNEG )  ! In/Out
           GAMMA(N) = BIG_GAMMA(N)

!if ( DO_EDDINGTON_FLUX ) then
!  write(*,'(A,i2,1p7e12.4)')' - ED Hom: gam1,gam2,biggam, eigenval,eigentrans,xneg(1:2)= ',&
!       N,Gamma1(n),GAMMA2(n), BIG_GAMMA(n), EIGENVALUE(n),EIGENTRANS(n), XNEG(1:2,N)
!else if ( DO_HEMISPHER_MEAN ) THEN
!  write(*,'(A,i2,1p7e12.4)')' - HM Hom: gam1,gam2,biggam, eigenval,eigentrans,xneg(1:2)=',&
!       N,Gamma1(n),GAMMA2(n), BIG_GAMMA(n), EIGENVALUE(n),EIGENTRANS(n), XNEG(1:2,N)
!endif

         else
           if ( CONSSCAT(N) ) THEN
             SAB(N)        = ZERO ; DAB(N)        = ZERO
             XPOS(1:2,N)   = ZERO ; XNEG(1:2,N)   = ZERO
             EIGENVALUE(N) = ZERO ; EIGENTRANS(n) = ZERO
             GAMMA(N)      = ZERO ; ALPHA(N)      = ZERO 
             EPSON(N)      = ZERO ; NORM_SAVED(N) = ZERO
           else
             CALL TWOSTREAM_HOM_SOLUTION &
             ( MAXLAYERS,                               & ! Dimensions
               N, FOURIER, STREAM_VALUE, PXSQ,          & ! Input
               OMEGA_TOTAL, ASYMM_TOTAL, DELTAU_VERT,   & ! Input
               SAB, DAB, EIGENVALUE, EIGENTRANS,        & ! In/Out
               XPOS, XNEG, ALPHA, EPSON, GAMMA, NORM_SAVED ) ! In/Out
           endif
         endif

!  Get Post-processing ("user") solutions for this layer
!   !@@ 2p3. 11/5/13. Post-processing control
!   4/12/24. Add consscat flag to each subroutine.
!   5/9/24. Will not be done for the FLux-only alternative methods (HemisphericMean, Eddington)

         IF ( DO_POSTPROCESSING ) THEN
           IF ( .not. DO_TRIDIAG ) THEN
             CALL TWOSTREAM_HOM_USERSOLUTION_PD &
              ( MAXLAYERS, MAX_USER_STREAMS, CONSSCAT(N),                & ! Dimensions/CSFlag
                N_USER_STREAMS, N, FOURIER, STREAM_VALUE, PX11,          & ! Input
                USER_STREAMS, ULP, XPOS, OMEGA_TOTAL, ASYMM_TOTAL,       & ! Input
                U_XPOS, U_XNEG, U_HELP_P, U_HELP_M )                       ! Output
           ELSE
             CALL TWOSTREAM_HOM_USERSOLUTION_TD &
              ( MAXLAYERS, MAX_USER_STREAMS, CONSSCAT(N),                 & ! Dimensions/CSFlag
                N_USER_STREAMS, N, FOURIER, STREAM_VALUE, PX11,           & ! Input
                USER_STREAMS, ULP, XPOS, OMEGA_TOTAL, ASYMM_TOTAL, GAMMA, & ! Input
                U_XPOS, U_XNEG, U_HELP_P, U_HELP_M )                        ! Output
           ENDIF
         ENDIF

!  end layer loop

      ENDDO

!  Prepare homogeneous solution multipliers
!   !@@ 2p3. 11/5/13. Post-processing control
!   !@@ 2p4. 8/15/14. User secants, Taylor-series control
!   5/9/24. Will not be done for the FLux-only alternative methods (HemisphericMean, Eddington)

      IF ( DO_POSTPROCESSING ) THEN
         CALL TWOSTREAM_HMULT_MASTER &
           ( MAXLAYERS, MAX_USER_STREAMS, CONSSCAT,                & ! Dimensions
             TAYLOR_ORDER, TAYLOR_SMALL, DELTAU_VERT, ASYMM_TOTAL, & ! Inputs 
             NLAYERS, N_USER_STREAMS, USER_SECANTS,  STREAM_VALUE, & ! Input
             EIGENVALUE, EIGENTRANS, T_DELT_USERM, TGMFUNC,        & ! Input
             ZETA_M, ZETA_P, HMULT_1, HMULT_2, CSPPHOMUP, CSPPHOMDN )   ! Output
      ENDIF

!  4.22.17. Boundary value problem - MATRIX PREPARATION (Tridiagonal solution)
!     Exception handling for Tridiagonal Matrix setup

!  check out # 12
!do n = 1, nlayers
!   write(*,*)n,EIGENVALUE(n),DELTAU_VERT(n),GAMMA(n),EIGENTRANS(n)
!enddo

      if ( DO_TRIDIAG ) then
        CALL TWOSTREAM_BVP_MATSETUP_TRIDIAG &
         ( MAXLAYERS, MAXTOTAL, DO_FULLQUAD,          & ! Input
           DO_INCLUDE_SURFACE, DO_BRDF_SURFACE,       & ! Input
           FOURIER, NLAYERS, NTOTAL, SURFACE_FACTOR,  & ! Input
           STREAM_VALUE, ALBEDO, BRDF_F, XPOS,        & ! Input
           GAMMA, EIGENTRANS, CONSSCAT, TGM,          & ! input
           AMAT_TD, BMAT_TD, DMAT_TD, EMAT_TD, SELM,  & ! output
           STATUS_SUB, MESSAGE )                        ! Output

!  Check out #14
!do n = 1, nlayers
! write(0,*) n,emat_td(n,1:4)
!enddo

        IF ( STATUS_SUB .NE. 0 ) THEN
          TRACE  = 'Call BVP_MATSETUP_TRIDIAG in 2S_FOURIER_MASTER'
          STATUS = 1 ; RETURN
        ENDIF
      endif

!  Boundary value problem - MATRIX PREPARATION (Pentadiagonal solution)
!     Pentadiagonal inverse option introduced, 25 June 2014. 4/22/17. Now controlled
!     Exception handling for Pentadiagonal Matrix setup

      if ( .not. DO_TRIDIAG ) then
        CALL TWOSTREAM_BVP_MATSETUP_PENTADIAG &
         ( MAXLAYERS, MAXTOTAL, DO_PENTADIAG_INVERSE, DO_FULLQUAD,  & ! Dimensions, flags
           DO_INCLUDE_SURFACE, DO_BRDF_SURFACE, FOURIER, NLAYERS,   & ! input
           NTOTAL, BVPSCALEFACTOR, SURFACE_FACTOR, STREAM_VALUE,    & ! input
           ALBEDO, BRDF_F, XPOS, XNEG, EIGENTRANS, CONSSCAT, TGM,   & ! input
           H_HOMP, H_HOMM, MAT, ELM, SELM,                          & ! output
           STATUS_SUB, MESSAGE )                                        ! Output
       IF ( STATUS_SUB .NE. 0 ) THEN
          TRACE  = 'Call BVP_MATSETUP_PENTADIAG in 2S_FOURIER_MASTER'
          STATUS = 1 ; RETURN
        ENDIF
      ENDIF

!  Thermal solutions
!     1. Find the Particular solution (NOT FOR transmittance only)
!     2. Compute thermal layer source terms. (Upwelling and Downwelling)
!       These will be scaled up by factor 4.pi if solar beams as well

!   !@@ 2p3. 11/5/13. Post-processing control
!   !@@ 2p4. 8/15/14. Greens function solution for Regular code
!   !@@ 2p4. 5/14/15. Thermal Cutoff variable introduced.
!  Rob Fix 5/24/17 for 2p4a, reintroduce classical PI subroutines
!           TWOSTREAM_THERMAL_CLSOLUTION
!           TWOSTREAM_THERMAL_CLSTERMS

!  4/12/24. Add Conservative Scattering Array, every subroutine
!  5/8/24. Introduce flag DO_HEMISPHER_MEAN for thermal regime
!    -- This has its own very simple solution. FLUXES ONLY, no postprocessing

      IF ( DO_INCLUDE_THERMEMISS ) THEN

        IF ( DO_HEMISPHER_MEAN ) THEN
          CALL TWOSTREAM_THERMAL_HMSOLUTION  &
              ( MAXLAYERS, NLAYERS, CONSSCAT, DELTAU_VERT, & ! Inputs
                THERMCOEFFS, GAMMA1, GAMMA2,               & ! Inputs
                T_WUPPER, T_WLOWER )                         ! Outputs
        ELSE
          IF ( DO_CLASSICAL_PI ) THEN
            CALL TWOSTREAM_THERMAL_CLSOLUTION  &
             ( MAXLAYERS, MAX_USER_STREAMS, DO_UPWELLING, DO_DNWELLING,       & ! Inputs
               CONSSCAT, NLAYERS, N_USER_STREAMS, STREAM_VALUE, USER_STREAMS, & ! Inputs
               OMEGA_TOTAL, ASYMM_TOTAL, SAB, DAB, DELTAU_VERT, THERMCOEFFS,  & ! Inputs
               T_WUPPER, T_WLOWER, U_TPOS2, U_TNEG2 )                           ! Outputs
          ELSE
            CALL TWOSTREAM_THERMAL_GFSOLUTION &
             ( MAXLAYERS, NLAYERS, OMEGA_TOTAL, DELTAU_VERT, THERMCOEFFS,   & ! Inputs
               CONSSCAT, TCUTOFF, EIGENVALUE, EIGENTRANS, XPOS, NORM_SAVED, & ! Input
               T_C_PLUS, T_C_MINUS, TTERM_SAVE, T_WUPPER, T_WLOWER )          ! Outputs
          ENDIF
        ENDIF

        IF ( DO_POSTPROCESSING.and..not. DO_HEMISPHER_MEAN) THEN
          IF ( DO_CLASSICAL_PI ) THEN
            CALL TWOSTREAM_THERMAL_CLSTERMS  &
              ( MAXLAYERS, MAX_USER_STREAMS, CONSSCAT,          & ! Inputs
                DO_UPWELLING, DO_DNWELLING, DO_SOLAR_SOURCES,   & ! Inputs
                NLAYERS, N_USER_STREAMS, PI4, USER_STREAMS,     & ! Inputs
                T_DELT_USERM, DELTAU_VERT, U_TPOS2, U_TNEG2,    & ! Inputs
                LAYER_TSUP_UP, LAYER_TSUP_DN  )                   ! Outputs
          ELSE
            CALL TWOSTREAM_THERMAL_GFSTERMS &
              ( MAXLAYERS, MAX_USER_STREAMS, CONSSCAT,          & ! Inputs
                DO_UPWELLING, DO_DNWELLING, DO_SOLAR_SOURCES,   & ! Input
                NLAYERS, N_USER_STREAMS, PI4, USER_STREAMS,     & ! Input
                TCUTOFF, T_DELT_USERM, DELTAU_VERT,             & ! Input
                U_XPOS, U_XNEG, HMULT_1, HMULT_2,               & ! Inputs
                T_C_PLUS, T_C_MINUS, TTERM_SAVE,                & ! Inputs
                LAYER_TSUP_UP, LAYER_TSUP_DN  )                   ! Output
          ENDIF
        ENDIF

      ENDIF

!write(*,*)eigenvalue(nlayers) ; pause
!      DO N = 1, NLAYERS
!         write(68,*)N,T_WLOWER(1:2,N)
!      ENDDO

!  Skip the thermal-only section if there are solar sources

      IF ( DO_SOLAR_SOURCES ) GO TO 455

!  ####################################################
!  Complete Radiation Field with Thermal-only solutions
!  ####################################################

!  Only one solution, local direct_beam flag NOT set

      IBEAM = 1
      DO_INCLUDE_DIRECTBEAM = .FALSE.

!  set the BVP PI solution at the lower/upper boundaries
!  2/12/24. No solutions if layer has Conservative scattering

      DO N = 1, NLAYERS
         DO I = 1, 2
            WUPPER(I,N) = T_WUPPER(I,N)
            WLOWER(I,N) = T_WLOWER(I,N)
         ENDDO
!         write(61,*)n,T_WLOWER(1:2,n),T_WUPPER(1:2,N)
!         write(62,*)n,WLOWER(1:2,n),T_WLOWER(1:2,N)
      ENDDO

!  Solve boundary value problem (Pentadiagonal solution). Add Control 4/22/17
!  Solve boundary value problem (Tridiagonal solution). 4/22/17

      IF ( DO_TRIDIAG ) THEN
        CALL TWOSTREAM_BVP_SOLUTION_TRIDIAG &
         ( MAXLAYERS, MAXBEAMS, MAXTOTAL, DO_FULLQUAD,        & ! Dimensions
           DO_INCLUDE_SURFACE, DO_INCLUDE_DIRECTBEAM,         & ! inputs
           DO_INCLUDE_SURFEMISS, DO_BRDF_SURFACE,             & ! inputs
           FOURIER, IBEAM, NLAYERS, NTOTAL, SURFACE_FACTOR,   & ! inputs
           STREAM_VALUE, ALBEDO, BRDF_F, EMISSIVITY, SURFBB,  & ! inputs
           DIRECT_BEAM, XPOS, XNEG, WUPPER, WLOWER, CONSSCAT, & ! inputs
           GAMMA, EIGENTRANS, AMAT_TD, BMAT_TD, DMAT_TD, EMAT_TD, SELM, & ! inputs
           LCON, MCON, LCON_XVEC, MCON_XVEC  )                  ! Output
      ENDIF

!     Version 2p3. Pentadiagonal inverse option introduced, 25 June 2014

      IF ( .not.DO_TRIDIAG ) THEN
        CALL TWOSTREAM_BVP_SOLUTION_PENTADIAG &
            ( MAXLAYERS, MAXBEAMS, MAXTOTAL,                   & ! Dimensions
              DO_PENTADIAG_INVERSE, DO_FULLQUAD,               & ! Flags
              DO_INCLUDE_SURFACE, DO_INCLUDE_DIRECTBEAM,       & ! inputs
              DO_INCLUDE_SURFEMISS, DO_BRDF_SURFACE,           & ! inputs
              FOURIER, IBEAM, NLAYERS, NTOTAL,                 & ! inputs
              BVPSCALEFACTOR, SURFACE_FACTOR, STREAM_VALUE,    & ! inputs
              ALBEDO, BRDF_F, EMISSIVITY, SURFBB,              & ! inputs
              DIRECT_BEAM, XPOS, XNEG, WUPPER, WLOWER,         & ! inputs
              MAT, ELM, SELM, CONSSCAT,                        & ! inputs
              H_PARTIC, LCON, MCON, LCON_XVEC, MCON_XVEC  )      ! Output
      ENDIF

!stop'BVPtest, thermal'

!  upwelling, MSMODE only, no Direct Beam inclusion.
!         !@@ 2p1 New OBSGEOM     option 12/21/12
!         !@@ 2p2 New 2S_LEVELOUT option 07/17/13
!         !@@ 2p3. 11/5/13. Post-processing control
!  Rob Fix 5/24/17 for 2p4a, reintroduce classical PI methods
!  Rob Fix 5/24/17 for 2p4a, Develop handles for Conservative scattering

      IF ( DO_UPWELLING .and. DO_POSTPROCESSING ) THEN
         CALL TWOSTREAM_UPUSER_INTENSITY &
           ( MAXLAYERS, MAXBEAMS, MAX_USER_STREAMS,                       & ! Dimensions
             DO_INCLUDE_SURFACE, DO_BRDF_SURFACE,                         & ! Surface control
             DO_CLASSICAL_PI, DO_TRIDIAG, DO_USER_OBSGEOMS,               & ! inputs !@@ 2p1, Control RT
             DO_SOLAR_SOURCES, DO_INCLUDE_THERMEMISS, DO_2S_LEVELOUT,     & ! inputs !@@ 2p2
             FOURIER, IBEAM, NLAYERS, N_USER_STREAMS, TAYLOR_ORDER,       & ! inputs !@@ 2p3 Greens
             BEAM_CUTOFF, PI4, SURFACE_FACTOR, ALBEDO, UBRDF_F,           & ! inputs
             FLUXMULT, STREAM_VALUE, TAYLOR_SMALL, DELTAU_VERT, CONSSCAT, & ! inputs
             GAMMA_P, GAMMA_M, SIGMA_P, ATERM_SAVE, BTERM_SAVE,           & ! Inputs !@@ 2p3 Greens
             INITIAL_TRANS, ITRANS_USERM, T_DELT_USERM, T_DELT_MUBAR,     & ! Inputs !@@ 2p3 Greens
             EIGENTRANS, LCON, LCON_XVEC, MCON, MCON_XVEC, WLOWER, TGM, CSPPHOMUP, & ! inputs
             U_XPOS, U_XNEG, U_WPOS2, HMULT_1, HMULT_2, EMULT_UP, LAYER_TSUP_UP,   & ! inputs
             IDOWNSURF, PMULT_UU, PMULT_UD,                               & ! Output !@@ 2p3 Greens
             INTENSITY_F_UP, RADLEVEL_F_UP, CUMSOURCE_UP )                  ! Output !@@ 2p2
      ENDIF

!  Downwelling, MSMODE only,
!         !@@ 2p1 New OBSGEOM     option 12/21/12
!         !@@ 2p2 New 2S_LEVELOUT option 07/17/13
!         !@@ 2p3. 11/5/13. Post-processing control
!  Rob Fix 5/24/17 for 2p4a, reintroduce classical PI methods
!  Rob Fix 5/24/17 for 2p4a, Develop handles for Conservative scattering

      IF ( DO_DNWELLING .and. DO_POSTPROCESSING ) THEN
         CALL TWOSTREAM_DNUSER_INTENSITY &
           ( MAXLAYERS, MAXBEAMS, MAX_USER_STREAMS,                      & ! Dimensions
             DO_CLASSICAL_PI, DO_TRIDIAG, DO_USER_OBSGEOMS,              & ! inputs !@@ 2p1, Control RT
             DO_SOLAR_SOURCES, DO_INCLUDE_THERMEMISS, DO_2S_LEVELOUT,    & ! inputs !@@ 2p2
             FOURIER, IBEAM, NLAYERS, N_USER_STREAMS, TAYLOR_ORDER, PI4, & ! inputs !@@ 2p3 Greens
             BEAM_CUTOFF, FLUXMULT, TAYLOR_SMALL, DELTAU_VERT, CONSSCAT, & ! inputs
             GAMMA_P, GAMMA_M, SIGMA_M, ATERM_SAVE, BTERM_SAVE,          & ! Inputs !@@ 2p3 Greens
             INITIAL_TRANS, ITRANS_USERM, T_DELT_USERM, T_DELT_MUBAR,    & ! Inputs !@@ 2p3 Greens
             LCON, MCON, LCON_XVEC, MCON_XVEC, U_XPOS, U_XNEG, U_WNEG2,  & ! Inputs
             HMULT_1, HMULT_2, EMULT_DN, LAYER_TSUP_DN, CSPPHOMDN,       & ! Inputs
             PMULT_DU, PMULT_DD,                                         & ! Output !@@ 2p3 Greens
             INTENSITY_F_DN, RADLEVEL_F_DN, CUMSOURCE_DN )                 ! Output !@@ 2p2
      ENDIF

!  Flux output. New Subroutine, 11/5/13 Version 2.3
!  Rob Fix 5/24/17 for 2p4a, Develop handles for Conservative scattering
!  Rob Fix 5/24/17 for 2p4a, Flux output at all levels
!  6/29/24. Include Eddington Flux flag in the argument list

      IF ( DO_INCLUDE_MVOUT ) THEN
         CALL TWOSTREAM_FLUXES &
            ( MAXBEAMS, MAXLAYERS, DO_UPWELLING, DO_DNWELLING,          & ! Dims/Flags
              DO_2S_LEVELOUT, DO_INCLUDE_DIRECTBEAM, DO_EDDINGTON_FLUX, & ! Flags
              IBEAM, NLAYERS, PI4, STREAM_VALUE, FLUX_FACTOR,           & ! Input Control
              FLUXMULT, X0, INITIAL_TRANS, TRANS_SOLAR_BEAM, CONSSCAT, TGM, & ! Input Control
              LCON, LCON_XVEC, MCON, MCON_XVEC, EIGENTRANS, WUPPER, WLOWER, & ! Input 2-stream solution
              FLUXES_TOA, FLUXES_BOA, FLUXLEVEL_UP, FLUXLEVEL_DN )            ! Output
      ENDIF

!  Finish Thermal only.

      RETURN

!  ##################################################
!  Complete Radiation Field with Solar Beam solutions
!  ##################################################

!  Continuation point

 455  CONTINUE

!  Start loop over various solar beams

      DO IBEAM = 1, NBEAMS

!  Solar beam Particular solution
!  ------------------------------

!  start layer loop

         DO N = 1, NLAYERS

!  Version 2p4 Greens function solution
!  Rob Fix 5/24/17 for 2p4a, reintroduce classical PI methods

!  5/9/24. Introduce alternative Particular integal for Eddington flux option
!    -- Possible to use average secant formulation in this.

            if ( DO_EDDINGTON_FLUX ) THEN
                CALL EDDINGTON_BEAM_SOLUTION &
                ( MAXLAYERS, MAXBEAMS, N, IBEAM, PI4, X0, FLUX_FACTOR,      & ! Inputs
                  BEAM_CUTOFF, AVERAGE_SECANT, INITIAL_TRANS, T_DELT_MUBAR, & ! Inputs
                  OMEGA_TOTAL, ASYMM_TOTAL, EIGENVALUE, GAMMA1, GAMMA2,     & ! Inputs
                  WVEC, WUPPER, WLOWER )                                      ! In/Out
            else
              IF ( DO_CLASSICAL_PI ) THEN
                CALL TWOSTREAM_CBEAM_SOLUTION &
                ( MAXLAYERS, MAXBEAMS, N, FOURIER, IBEAM, PI4,    & ! Inputs
                  FLUX_FACTOR, BEAM_CUTOFF, STREAM_VALUE, PX0X,   & ! Inputs
                  AVERAGE_SECANT, INITIAL_TRANS, T_DELT_MUBAR,    & ! Inputs
                  OMEGA_TOTAL, ASYMM_TOTAL, SAB, DAB, EIGENVALUE, & ! Inputs
                  QSUMVEC, QDIFVEC, QVEC, WVEC, WUPPER, WLOWER )    ! In/Out
              ELSE
                CALL TWOSTREAM_GBEAM_SOLUTION &
                ( MAXLAYERS, MAXBEAMS,                                & ! Dimensions
                  TAYLOR_ORDER, TAYLOR_SMALL, DELTAU_VERT,            & ! Inputs 
                  N, FOURIER, IBEAM, PI4, FLUX_FACTOR,                & ! Inputs
                  BEAM_CUTOFF, PX0X, OMEGA_TOTAL, ASYMM_TOTAL,        & ! Inputs
                  AVERAGE_SECANT, INITIAL_TRANS, T_DELT_MUBAR,        & ! Inputs
                  XPOS, EIGENVALUE, EIGENTRANS, NORM_SAVED,           & ! Inputs
                  GAMMA_M, GAMMA_P, DMI, DPI, ATERM_SAVE, BTERM_SAVE, & ! Output
                  CFUNC, DFUNC, GFUNC_UP, GFUNC_DN, WUPPER, WLOWER )    ! Output
              ENDIF
            endif

!  Version 2p4 Classical solution
!  5/9/24. Not for Eddington flux option

            if ( .not. DO_EDDINGTON_FLUX ) then
              if ( CONSSCAT(N) ) then
                CALL TWOSTREAM_CONSSCAT_SOLUTION &
                ( MAXLAYERS, MAXBEAMS, DO_PLANE_PARALLEL, N, IBEAM, PI4, X0(IBEAM), & ! Inputs
                  FLUX_FACTOR, BEAM_CUTOFF, STREAM_VALUE, ASYMM_TOTAL,              & ! Inputs
                  INITIAL_TRANS, AVERAGE_SECANT, T_DELT_MUBAR,                      & ! Inputs
                  WVEC, WUPPER, WLOWER )                                              ! In/Out
              endif
            endif

!  FOLLOWING CODE DOES NOT APPLY for GREENS FUNCTION (VERSION 2p4)
!    Classical-PI User solutions. Version 2p4a, re-introduced
!   5/9/24. Will not be done for the Flux-only, Eddington method

            IF (  DO_POSTPROCESSING.and.DO_CLASSICAL_PI ) THEN
               CALL TWOSTREAM_CBEAM_USERSOLUTION &
                ( DO_UPWELLING, DO_DNWELLING,                         & ! Input
                  NLAYERS, NBEAMS, N_USER_STREAMS, N, FOURIER, IBEAM, & ! Input
                  FLUX_FACTOR, BEAM_CUTOFF, STREAM_VALUE, PI4, PX11,  & ! Input
                  OMEGA_TOTAL, ASYMM_TOTAL, USER_STREAMS, ULP, WVEC,  & ! Input
                  U_WPOS2, U_WNEG2, W_HELP )                            ! Output
            ENDIF

!  end layer loop

         END DO

!  Add thermal solutions if flagged. NO modulus on the thermal contribution.

         IF ( DO_INCLUDE_THERMEMISS ) THEN
            DO N = 1, NLAYERS
               DO I = 1, 2
                  WUPPER(I,N) = WUPPER(I,N) + T_WUPPER(I,N)
                  WLOWER(I,N) = WLOWER(I,N) + T_WLOWER(I,N)
               ENDDO
            ENDDO
         ENDIF

!  Solve boundary value problem (Tridiagonal solution). 4/22/17

         IF ( DO_TRIDIAG ) THEN
           CALL TWOSTREAM_BVP_SOLUTION_TRIDIAG &
            ( MAXLAYERS, MAXBEAMS, MAXTOTAL, DO_FULLQUAD,        & ! Dimensions
              DO_INCLUDE_SURFACE, DO_DIRECTBEAM(IBEAM),          & ! inputs
              DO_INCLUDE_SURFEMISS, DO_BRDF_SURFACE,             & ! inputs
              FOURIER, IBEAM, NLAYERS, NTOTAL, SURFACE_FACTOR,   & ! inputs
              STREAM_VALUE, ALBEDO, BRDF_F, EMISSIVITY, SURFBB,  & ! inputs
              DIRECT_BEAM, XPOS, XNEG, WUPPER, WLOWER, CONSSCAT, & ! inputs
              GAMMA, EIGENTRANS, AMAT_TD, BMAT_TD, DMAT_TD, EMAT_TD, SELM, & ! inputs
              LCON, MCON, LCON_XVEC, MCON_XVEC  )                  ! Output
         ENDIF

!  Solve boundary value problem (Pentadiagonal solution). 4/22/17 Control added
!     Pentadiagonal inverse option introduced, 25 June 2014

         IF ( .not.DO_TRIDIAG ) THEN
           CALL TWOSTREAM_BVP_SOLUTION_PENTADIAG &
            ( MAXLAYERS, MAXBEAMS, MAXTOTAL,                   & ! Dimensions
              DO_PENTADIAG_INVERSE, DO_FULLQUAD,               & ! Flags
              DO_INCLUDE_SURFACE, DO_DIRECTBEAM(IBEAM),        & ! inputs
              DO_INCLUDE_SURFEMISS, DO_BRDF_SURFACE,           & ! inputs
              FOURIER, IBEAM, NLAYERS, NTOTAL,                 & ! inputs
              BVPSCALEFACTOR, SURFACE_FACTOR, STREAM_VALUE,    & ! inputs
              ALBEDO, BRDF_F, EMISSIVITY, SURFBB,              & ! inputs
              DIRECT_BEAM, XPOS, XNEG, WUPPER, WLOWER,         & ! inputs
              MAT, ELM, SELM, CONSSCAT,                        & ! inputs
              H_PARTIC, LCON, MCON, LCON_XVEC, MCON_XVEC  )      ! Output
         ENDIF

!stop'BVPtest, Solar'

! ##################################
!   Radiance Field Post Processing
! ##################################

!  upwelling, MSMODE only, no Direct Beam inclusion.
!         !@@ 2p1 New OBSGEOM     option 12/21/12
!         !@@ 2p2 New 2S_LEVELOUT option 07/17/13
!         !@@ 2p3. 11/5/13. Post-processing control
!  Rob Fix 5/24/17 for 2p4a, reintroduce classical PI methods
!  Rob Fix 5/24/17 for 2p4a, Develop handles for Conservative scattering

!  Rob Fix. This flag should be set true if you want to include the Direct Downwelling Flux.

         DO_INCLUDE_DIRECTBEAM = .TRUE.

         IF ( DO_UPWELLING .and. DO_POSTPROCESSING ) THEN
            CALL TWOSTREAM_UPUSER_INTENSITY &
              ( MAXLAYERS, MAXBEAMS, MAX_USER_STREAMS,                       & ! Dimensions
                DO_INCLUDE_SURFACE, DO_BRDF_SURFACE,                         & ! Surface control
                DO_CLASSICAL_PI, DO_TRIDIAG, DO_USER_OBSGEOMS,               & ! inputs !@@ 2p1, Control RT
                DO_SOLAR_SOURCES, DO_INCLUDE_THERMEMISS, DO_2S_LEVELOUT,     & ! inputs !@@ 2p2
                FOURIER, IBEAM, NLAYERS, N_USER_STREAMS, TAYLOR_ORDER,       & ! inputs !@@ 2p3 Greens
                BEAM_CUTOFF, PI4, SURFACE_FACTOR, ALBEDO, UBRDF_F,           & ! inputs
                FLUXMULT, STREAM_VALUE, TAYLOR_SMALL, DELTAU_VERT, CONSSCAT, & ! inputs
                GAMMA_P, GAMMA_M, SIGMA_P, ATERM_SAVE, BTERM_SAVE,           & ! Inputs !@@ 2p3 Greens
                INITIAL_TRANS, ITRANS_USERM, T_DELT_USERM, T_DELT_MUBAR,     & ! Inputs !@@ 2p3 Greens
                EIGENTRANS, LCON, LCON_XVEC, MCON, MCON_XVEC, WLOWER, TGM, CSPPHOMUP, & ! inputs
                U_XPOS, U_XNEG, U_WPOS2, HMULT_1, HMULT_2, EMULT_UP, LAYER_TSUP_UP,   & ! inputs
                IDOWNSURF, PMULT_UU, PMULT_UD,                               & ! Output !@@ 2p3 Greens
                INTENSITY_F_UP, RADLEVEL_F_UP, CUMSOURCE_UP )                  ! Output !@@ 2p2
        ENDIF

!  Downwelling, MSMODE only,
!         !@@ 2p1 New OBSGEOM     option 12/21/12
!         !@@ 2p2 New 2S_LEVELOUT option 07/17/13
!         !@@ 2p3. 11/5/13. Post-processing control
!  Rob Fix 5/24/17 for 2p4a, reintroduce classical PI methods
!  Rob Fix 5/24/17 for 2p4a, Develop handles for Conservative scattering

         IF ( DO_DNWELLING .and. DO_POSTPROCESSING ) THEN
            CALL TWOSTREAM_DNUSER_INTENSITY &
              ( MAXLAYERS, MAXBEAMS, MAX_USER_STREAMS,                      & ! Dimensions
                DO_CLASSICAL_PI, DO_TRIDIAG, DO_USER_OBSGEOMS,              & ! inputs !@@ 2p1, Control RT
                DO_SOLAR_SOURCES, DO_INCLUDE_THERMEMISS, DO_2S_LEVELOUT,    & ! inputs !@@ 2p2
                FOURIER, IBEAM, NLAYERS, N_USER_STREAMS, TAYLOR_ORDER, PI4, & ! inputs !@@ 2p3 Greens
                BEAM_CUTOFF, FLUXMULT, TAYLOR_SMALL, DELTAU_VERT, CONSSCAT, & ! inputs
                GAMMA_P, GAMMA_M, SIGMA_M, ATERM_SAVE, BTERM_SAVE,          & ! Inputs !@@ 2p3 Greens
                INITIAL_TRANS, ITRANS_USERM, T_DELT_USERM, T_DELT_MUBAR,    & ! Inputs !@@ 2p3 Greens
                LCON, MCON, LCON_XVEC, MCON_XVEC, U_XPOS, U_XNEG, U_WNEG2,  & ! Inputs
                HMULT_1, HMULT_2, EMULT_DN, LAYER_TSUP_DN, CSPPHOMDN,       & ! Inputs
                PMULT_DU, PMULT_DD,                                         & ! Output !@@ 2p3 Greens
                INTENSITY_F_DN, RADLEVEL_F_DN, CUMSOURCE_DN )                 ! Output !@@ 2p2
         ENDIF

!  Flux output. New Subroutine, 11/5/13 Version 2.3
!  Rob Fix 5/24/17 for 2p4a, Develop handles for Conservative scattering
!  Rob Fix 5/24/17 for 2p4a, Flux output at all levels
!  6/29/24. Include Eddington Flux flag in the argument list

         IF ( DO_INCLUDE_MVOUT ) THEN
            CALL TWOSTREAM_FLUXES &
            ( MAXBEAMS, MAXLAYERS, DO_UPWELLING, DO_DNWELLING,          & ! Dims/Flags
              DO_2S_LEVELOUT, DO_INCLUDE_DIRECTBEAM, DO_EDDINGTON_FLUX, & ! Flags
              IBEAM, NLAYERS, PI4, STREAM_VALUE, FLUX_FACTOR,           & ! Input Control
              FLUXMULT, X0, INITIAL_TRANS, TRANS_SOLAR_BEAM, CONSSCAT, TGM, & ! Input Control
              LCON, LCON_XVEC, MCON, MCON_XVEC, EIGENTRANS, WUPPER, WLOWER, & ! Input 2-stream solution
              FLUXES_TOA, FLUXES_BOA, FLUXLEVEL_UP, FLUXLEVEL_DN )            ! Output
         ENDIF

!  End loop over beam solutions

      END DO

!  ######
!  finish
!  ######

      RETURN
END SUBROUTINE TWOSTREAM_FOURIER_MASTER

end module twostream_master_m
