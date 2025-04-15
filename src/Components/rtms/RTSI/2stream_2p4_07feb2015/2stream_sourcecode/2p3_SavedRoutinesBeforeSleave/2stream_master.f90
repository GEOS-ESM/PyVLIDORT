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

Use twostream_miscsetups_m
Use twostream_solutions_m
Use twostream_bvproblem_m
Use twostream_intensity_m
Use twostream_thermalsup_m

PUBLIC

contains

SUBROUTINE TWOSTREAM_MASTER &
        ( MAXLAYERS, MAXTOTAL, MAXTHREADS, MAXMESSAGES,                   & ! Dimensions
          MAXBEAMS, MAX_USER_STREAMS, MAX_USER_RELAZMS,                   & ! Dimensions
          MAX_GEOMETRIES, MAX_USER_OBSGEOMS,                              & ! Dimensions !@@ 2p1
          DO_UPWELLING, DO_DNWELLING, DO_PLANE_PARALLEL, DO_2S_LEVELOUT,  & ! Inputs     !@@ 2p2
          DO_MVOUT_ONLY, DO_ADDITIONAL_MVOUT,                             & ! Inputs     !@@ 2p3
          DO_SOLAR_SOURCES, DO_THERMAL_EMISSION, DO_SURFACE_EMISSION,     & ! Inputs
          DO_D2S_SCALING, DO_BRDF_SURFACE, DO_USER_OBSGEOMS,              & ! Inputs     !@@ 2p1
          THREAD, NLAYERS, NTOTAL, STREAM_VALUE,                          & ! Inputs
          N_USER_OBSGEOMS, USER_OBSGEOMS,                                 & ! Inputs     !@@ 2p1
          N_USER_STREAMS, USER_ANGLES, N_USER_RELAZMS, USER_RELAZMS,      & ! Inputs
          FLUX_FACTOR, NBEAMS, BEAM_SZAS, EARTH_RADIUS, HEIGHT_GRID,      & ! Inputs
          DELTAU_INPUT, OMEGA_INPUT, ASYMM_INPUT, D2S_SCALING,            & ! Inputs
          THERMAL_BB_INPUT, LAMBERTIAN_ALBEDO, BRDF_F_0, BRDF_F, UBRDF_F, & ! Inputs
          EMISSIVITY, SURFBB,                                             & ! Inputs
          INTENSITY_TOA, INTENSITY_BOA, FLUXES_TOA, FLUXES_BOA,           & ! Outputs !@@ 2p3(Fluxes)
          RADLEVEL_UP, RADLEVEL_DN, N_GEOMETRIES,                         & ! Outputs !@@ 2p2
          STATUS_INPUTCHECK, C_NMESSAGES, C_MESSAGES, C_ACTIONS,          & ! Exception handling
          STATUS_EXECUTION,  E_MESSAGE, E_TRACE_1, E_TRACE_2 )              ! Exception handling

      implicit none

!  Precision

      INTEGER, PARAMETER :: dp     = KIND( 1.0D0 )

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

!  Subroutine input arguments
!  --------------------------

!  Dimensions :
!      MAXTOTAL       = 2 * MAXLAYERS
!      MAX_GEOMETRIES = MAXBEAMS * MAX_USER_STREAMS * MAX_USER_RELAZMS
!      !@@ MAX_USER_OBSGEOMS >/= MAXBEAMS

      INTEGER, INTENT(IN)        :: MAXTHREADS, MAXMESSAGES
      INTEGER, INTENT(IN)        :: MAXLAYERS, MAXTOTAL
      INTEGER, INTENT(IN)        :: MAXBEAMS, MAX_GEOMETRIES, MAX_USER_OBSGEOMS  !@@ 2p1
      INTEGER, INTENT(IN)        :: MAX_USER_STREAMS, MAX_USER_RELAZMS

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

!  Input thread

      INTEGER, INTENT(IN)        :: THREAD

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

      REAL(kind=dp), INTENT(IN)  :: FLUX_FACTOR

!  Solar geometry. [Now Intent(inout), thanks to option for ObsGeom !@@ 2p1

      INTEGER, INTENT(INOUT)        :: NBEAMS
      REAL(kind=dp), INTENT(INOUT)  :: BEAM_SZAS ( MAXBEAMS )

!  Height and earth radius (latter could be re-set internally)

      REAL(kind=dp), INTENT(INOUT) :: EARTH_RADIUS
      REAL(kind=dp), INTENT(IN)    :: HEIGHT_GRID ( 0:MAXLAYERS )

!  Geometry specification height
!      REAL(kind=dp), INTENT(IN)  :: GEOMETRY_SPECHEIGHT

!  Atmospheric optical properties

      REAL(kind=dp), INTENT(IN)  :: DELTAU_INPUT(MAXLAYERS, MAXTHREADS)
      REAL(kind=dp), INTENT(IN)  :: OMEGA_INPUT (MAXLAYERS, MAXTHREADS)
      REAL(kind=dp), INTENT(IN)  :: ASYMM_INPUT (MAXLAYERS, MAXTHREADS)
      REAL(kind=dp), INTENT(IN)  :: D2S_SCALING (MAXLAYERS, MAXTHREADS)

!  Atmospheric thermal sources

      REAL(kind=dp), INTENT(IN)  :: THERMAL_BB_INPUT ( 0:MAXLAYERS )

!  Lambertian surface control (threaded)

      REAL(kind=dp), INTENT(IN)  :: LAMBERTIAN_ALBEDO (MAXTHREADS)

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
      REAL(kind=dp), INTENT(IN)  :: SURFBB

!  Output
!  ------

!  Radiance Results

      REAL(kind=dp), INTENT(INOUT) :: INTENSITY_TOA(MAX_GEOMETRIES,MAXTHREADS)
      REAL(kind=dp), INTENT(INOUT) :: INTENSITY_BOA(MAX_GEOMETRIES,MAXTHREADS)

!  Flux output
!     ! @@ Rob Spurr, 05 November 2013, Version 2.3 --> Flux Output

     REAL(kind=dp), INTENT(INOUT) :: FLUXES_TOA(MAXBEAMS,2,MAXTHREADS)
     REAL(kind=dp), INTENT(INOUT) :: FLUXES_BOA(MAXBEAMS,2,MAXTHREADS)

!  output solutions at ALL levels
!     ! @@ Rob Spurr, 17 July 2013, Version 2.2 --> Optional Output at ALL LEVELS

      REAL(kind=dp), INTENT(INOUT) :: RADLEVEL_UP (MAX_GEOMETRIES,0:MAXLAYERS,MAXTHREADS)
      REAL(kind=dp), INTENT(INOUT) :: RADLEVEL_DN (MAX_GEOMETRIES,0:MAXLAYERS,MAXTHREADS)

!  Numbers (geometry)
!   N_GEOMETRIES = NBEAMS * N_USER_STREAMS * N_USER_RELAZMS (Lattice value)
!   N_GEOMETRIES = N_USER_OBSGEOMS                          (OBsGeom value)

      INTEGER, INTENT(INOUT)       :: N_GEOMETRIES

!  Exception handling
!  ------------------

!    1. Check Messages and actions

      INTEGER      , INTENT(OUT) :: STATUS_INPUTCHECK
      INTEGER      , INTENT(OUT) :: C_NMESSAGES
      CHARACTER*100, INTENT(OUT) :: C_MESSAGES(MAXMESSAGES)
      CHARACTER*100, INTENT(OUT) :: C_ACTIONS (MAXMESSAGES)

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
      REAL(kind=dp) :: ASYMM_TOTAL(MAXLAYERS)

!  Chapman factors (from pseudo-spherical geometry)

      REAL(kind=dp) :: CHAPMAN_FACTORS ( MAXLAYERS, MAXLAYERS, MAXBEAMS )
      REAL(kind=dp) :: LOCAL_SZA       ( 0:MAXLAYERS, MAXBEAMS )

!     Last layer to include Particular integral solution
!     Average-secant and initial tramsittance factors for solar beams.
!     Solar beam attenuation

      INTEGER       :: LAYER_PIS_CUTOFF ( MAXBEAMS )
      REAL(kind=dp) :: INITIAL_TRANS    ( MAXLAYERS, MAXBEAMS )
      REAL(kind=dp) :: AVERAGE_SECANT   ( MAXLAYERS, MAXBEAMS )
      REAL(kind=dp) :: LOCAL_CSZA       ( MAXLAYERS, MAXBEAMS )
      REAL(kind=dp) :: SOLAR_BEAM_OPDEP ( MAXBEAMS )

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

!  Forcing term multipliers (saved for whole atmosphere)

      REAL(kind=dp) :: EMULT_UP (MAX_USER_STREAMS,MAXLAYERS,MAXBEAMS)
      REAL(kind=dp) :: EMULT_DN (MAX_USER_STREAMS,MAXLAYERS,MAXBEAMS)

!  Fourier-component solutions

      REAL(kind=dp) :: INTENSITY_F_UP (MAX_USER_STREAMS,MAXBEAMS)
      REAL(kind=dp) :: INTENSITY_F_DN (MAX_USER_STREAMS,MAXBEAMS)

!  Fourier-component solutions at ALL levels
!     ! @@ Rob Spurr, 17 July 2013, Version 2.2 --> Optional Output at ALL LEVELS

      REAL(kind=dp) :: RADLEVEL_F_UP (MAX_USER_STREAMS,MAXBEAMS,0:MAXLAYERS)
      REAL(kind=dp) :: RADLEVEL_F_DN (MAX_USER_STREAMS,MAXBEAMS,0:MAXLAYERS)

!  Single scatter solutions, commented out in this  version
!      REAL(kind=dp) :: INTENSITY_SS_UP(N_GEOMETRIES)
!      REAL(kind=dp) :: INTENSITY_SS_DN(N_GEOMETRIES)

!  Other local variables
!  ---------------------

!  Local error handling

      CHARACTER(LEN=3) :: CF, WTHREAD
      LOGICAL          :: DO_USER_STREAMS
      INTEGER          :: FOURIER, N_FOURIERS, STATUS_SUB
      INTEGER          :: N, UA, UM, IB, N_VIEWING, IBEAM, I, T, LUM, LUA
      REAL(kind=dp)    :: AZM_ARGUMENT, DFC, DEG_TO_RAD, PI4
      REAL(kind=dp)    :: OMFAC, M1FAC, GDIFF, ALBEDO

!  Geometry offset arrays

      INTEGER          :: IBOFF ( MAXBEAMS )
      INTEGER          :: UMOFF ( MAXBEAMS, MAX_USER_STREAMS )

!  Local azimuth factors

      REAL(kind=dp)    :: AZMFAC (MAX_USER_STREAMS,MAXBEAMS,MAX_USER_RELAZMS)

!  Post-processing flag (new for Version 2p3)

      LOGICAL          :: DO_POSTPROCESSING

!  Cosines and sines

      REAL(kind=dp)    :: X0  ( MAXBEAMS )
      REAL(kind=dp)    :: USER_STREAMS ( MAX_USER_STREAMS )
      REAL(kind=dp)    :: MUSTREAM, SINSTREAM

!  Thermal help variables

      REAL(kind=dp)    :: TCOM1 ( MAXLAYERS, 2 )
      REAL(kind=dp)    :: DELTAU_POWER ( MAXLAYERS, 2 )

!mick - singularity buster output
      LOGICAL          :: SBUST(6)

!  Thread number

      WTHREAD = '000'
      IF (THREAD.LT.10)WRITE(WTHREAD(3:3),'(I1)')THREAD
      IF (THREAD.GT.99)WRITE(WTHREAD(1:3),'(I3)')THREAD
      IF (THREAD.GE.10.and.THREAD.LE.99)WRITE(WTHREAD(2:3),'(I2)')THREAD
      T = THREAD

!  Local user indices; !@@ Only required for OBSGEOM option

      LUM = 1
      LUA = 1

!  Initialize output arrays for current thread
!  -------------------------------------------

!mick fix 11/8/2012 - added main output
!  Main output

      INTENSITY_TOA(:,T) = 0.0d0
      INTENSITY_BOA(:,T) = 0.0d0

! @@ Rob Spurr, 17 July 2013, Version 2.2 --> Optional Output at ALL LEVELS

      RADLEVEL_UP (:,:,T) = 0.0d0
      RADLEVEL_DN (:,:,T) = 0.0d0

! @@ Rob Spurr, 05 November 2013, Version 2.3 --> BOA_TOA Flux outputs

      FLUXES_TOA(:,:,T) = 0.0d0
      FLUXES_BOA(:,:,T) = 0.0d0

!  Input check

      STATUS_INPUTCHECK = 0
      C_NMESSAGES       = 0

!  Execution status and message/traces

      STATUS_EXECUTION  = 0
      E_MESSAGE = ' '
      E_TRACE_1 = ' '
      E_TRACE_2 = ' '

!  Constants
!  ---------

      DEG_TO_RAD = DACOS(-1.0d0)/180.0d0
      PI4 = DEG_TO_RAD * 720.0d0

!  Input checking
!  ==============

!  Check input Basic. This could be put outside the thread loop.
!    SS inputs are omitted in this version........
!    !@@ 2p1, Observational Geometry inputs are included (New 12/21/12)
!    !@@ 2p3 Includes check on Flux output flags, and setting of Post-Processing flag

      CALL TWOSTREAM_CHECK_INPUTS_BASIC  &
        ( MAXLAYERS, MAXMESSAGES, MAX_USER_OBSGEOMS,             & ! Dimensions !@@
          MAXBEAMS, MAX_USER_STREAMS, MAX_USER_RELAZMS,          & ! Dimensions
          DO_UPWELLING, DO_DNWELLING, DO_PLANE_PARALLEL,         & ! Input
          DO_SOLAR_SOURCES, DO_THERMAL_EMISSION,                 & ! Input
          DO_MVOUT_ONLY, DO_ADDITIONAL_MVOUT, DO_POSTPROCESSING, & ! Input !@@ New line, 2p3
          DO_USER_OBSGEOMS, N_USER_OBSGEOMS, USER_OBSGEOMS,      & ! Input !@@ New line
          NLAYERS, NBEAMS, N_USER_STREAMS, N_USER_RELAZMS,       & ! Input
          BEAM_SZAS, USER_ANGLES, USER_RELAZMS,                  & ! Input
          EARTH_RADIUS, HEIGHT_GRID,                             & ! Input
          STATUS_SUB, C_NMESSAGES, C_MESSAGES, C_ACTIONS )         ! Output

      IF ( STATUS_SUB .EQ. 1 ) THEN
        STATUS_INPUTCHECK = 1
        RETURN
      ENDIF

!  Check input threaded values (IOPs in the atmosphere)

     CALL TWOSTREAM_CHECK_INPUTS_THREAD &
       ( MAXLAYERS, MAXTHREADS, MAXMESSAGES,               & ! Dimensions
         NLAYERS, THREAD,                                  & ! input
         DELTAU_INPUT, OMEGA_INPUT, ASYMM_INPUT,           & ! Input
         STATUS_SUB, C_NMESSAGES, C_MESSAGES, C_ACTIONS )    ! Output

      IF ( STATUS_SUB .EQ. 1 ) THEN
        STATUS_INPUTCHECK = 1
        RETURN
      ENDIF

!  Geometry offsets
!  ================

!  Save some offsets for indexing geometries

!   !@@ 2p1, This section revised for the Observational Geometry option
!   !@@ N_GEOMETRIES = NBEAMS * N_USER_STREAMS * N_USER_RELAZMS
!   !@@ 2p3, This section revised for the post-processing flag

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

!  Quadrature

      MUSTREAM  = STREAM_VALUE
      SINSTREAM = DSQRT(1.0d0-MUSTREAM*MUSTREAM)

!  Solar zenith angle cosine

      DO IB = 1, NBEAMS
        X0(IB) = DCOS ( BEAM_SZAS(IB) * DEG_TO_RAD )
      ENDDO

!  User stream cosines. 11/5/13 2p3 Post-processing control

      IF ( DO_POSTPROCESSING ) THEN
         DO I = 1, N_USER_STREAMS
            USER_STREAMS(I) = DCOS(DEG_TO_RAD * USER_ANGLES(I))
         ENDDO
      ENDIF

!  Set local atmospheric optical properties (Apply delta 2s scaling)
!  Just copy inputs, if not required

      IF ( DO_D2S_SCALING ) THEN
        DO N = 1, NLAYERS
          !OMFAC = 1.0D0 - OMEGA(I)*FA (= alpha)
          OMFAC = 1.0d0 - OMEGA_INPUT(N,T) * D2S_SCALING(N,T)
          !M1FAC = 1.0D0 - FA
          M1FAC = 1.0d0 - D2S_SCALING(N,T)
          !GDIFF = PF(I,1) - FA
          GDIFF = ASYMM_INPUT(N,T) - D2S_SCALING(N,T)
          DELTAU_VERT(N) = OMFAC * DELTAU_INPUT(N,T)
          OMEGA_TOTAL(N) = M1FAC * OMEGA_INPUT(N,T) / OMFAC
          ASYMM_TOTAL(N) = GDIFF / M1FAC
        ENDDO
      ELSE
        DO N = 1, NLAYERS
          DELTAU_VERT(N) = DELTAU_INPUT(N,T)
          OMEGA_TOTAL(N) = OMEGA_INPUT(N,T)
          ASYMM_TOTAL(N) = ASYMM_INPUT(N,T)
        ENDDO
      ENDIF

!mick fix 1/7/2012 - singularity busters added

!  Note: Due to delta-m scaling, omega and/or g may be
!        modified in such a way as to make them unphysical or introduce
!        instability in the two-stream case; therefore, we recheck omega
!        and g AFTER delta-m scaling and slightly adjust them if necessary

      DO N = 1, NLAYERS
        SBUST = .false.

        !Singularity buster for single scatter albedo
        IF (OMEGA_TOTAL(N) > 0.999999999D0) THEN
          OMEGA_TOTAL(N) = 0.999999999D0
          SBUST(1) = .true.
        ELSE IF (OMEGA_TOTAL(N) < 1.0D-9) THEN
          OMEGA_TOTAL(N) = 1.0D-9
          SBUST(2) = .true.
        END IF

        !Singularity buster for asymmetry parameter
        IF (ASYMM_TOTAL(N) > 0.999999999D0) THEN
          ASYMM_TOTAL(N) = 0.999999999D0
          SBUST(3) = .true.
        ELSE IF (ASYMM_TOTAL(N) < -0.999999999D0) THEN
          ASYMM_TOTAL(N) = -0.999999999D0
          SBUST(4) = .true.
        ELSE IF ((ASYMM_TOTAL(N) >= 0.0D0) .AND. &
                 (ASYMM_TOTAL(N) < 1.0D-9)) THEN
          ASYMM_TOTAL(N) = 1.0D-9
          SBUST(5) = .true.
        ELSE IF ((ASYMM_TOTAL(N) < 0.0D0) .AND. &
                 (ASYMM_TOTAL(N) > -1.0D-9)) THEN
          ASYMM_TOTAL(N) = -1.0D-9
          SBUST(6) = .true.
        END IF

        !WRITE(*,*)
        !WRITE(*,'(A,I2)') 'FOR LAYER: ',N
        !DO I=1,6
        !  WRITE(*,'(A,I1,A,L1)') '  SBUST(',I,') = ',SBUST(I)
        !ENDDO

      ENDDO

!  Initialise Fourier loop
!  =======================

!  Set Fourier number, Nominally 1 in absence of SS-only flag
!    Zero if no solar sources (Thermal-only run)
!  !@@ 2p3, Set NFOURIERS equal to zero for MVOUT_ONLY

      N_FOURIERS = 1
      IF (  DO_MVOUT_ONLY )         N_FOURIERS = 0
      IF ( .NOT. DO_SOLAR_SOURCES ) N_FOURIERS = 0

!mick fix 1/7/2012 - (test - make this permanent?)
      IF ( (NBEAMS == 1) .AND. (BEAM_SZAS(1) < 1.0D-8) ) &
        N_FOURIERS = 0

!  Albedo

      ALBEDO = LAMBERTIAN_ALBEDO(T)

!  Old code was dependent on SS flag
!      IF ( DO_SSFULL ) THEN
!        N_FOURIERS = 0
!      ELSE
!        N_FOURIERS = 1
!      ENDIF

!  Fourier loop
!  ============

      DO FOURIER = 0, N_FOURIERS

!  Azimuth cosine factors (Fourier = 1). !@@ 2p1, Notice OBSGEOM option
!  !@@ 2p3, not required for FLux-only output

        AZMFAC = 0.0D0
        IF ( DO_POSTPROCESSING ) THEN
          IF ( FOURIER .GT. 0 ) THEN
            DFC = DBLE(FOURIER)
            IF ( DO_USER_OBSGEOMS.and.DO_SOLAR_SOURCES ) THEN
              DO IB = 1, NBEAMS
                AZM_ARGUMENT = USER_RELAZMS(IB) * DFC
                AZMFAC(LUM,IB,LUA) = DCOS(DEG_TO_RAD*AZM_ARGUMENT)
              ENDDO
            ELSE
              DO UA = 1, N_USER_RELAZMS
                DO IB = 1, NBEAMS
                  DO UM = 1, N_USER_STREAMS
                    AZM_ARGUMENT = USER_RELAZMS(UA) * DFC
                    AZMFAC(UM,IB,UA) = DCOS(DEG_TO_RAD*AZM_ARGUMENT)
                  ENDDO
                ENDDO
              ENDDO
            ENDIF
          ENDIF
        ENDIF

!  Main call to Lidort Fourier module.
!  ----------------------------------

!  !@@ Add Observational Geometry dimension and control variables  !@@ 2p1
!  !@@ Call statement expanded to include ALL-LEVEL outputs        !@@ 2p2
!  !@@ Call statement expanded to include Flux outputs and control !@@ 2p3  11/5/13

        CALL TWOSTREAM_FOURIER_MASTER &
        ( MAXLAYERS, MAXTOTAL, MAXBEAMS, MAX_USER_STREAMS, MAXTHREADS,   & ! Dimensions
          DO_UPWELLING, DO_DNWELLING, DO_BRDF_SURFACE, DO_USER_OBSGEOMS, & ! Input flags control (including Obsgeom)
          DO_SOLAR_SOURCES, DO_THERMAL_EMISSION, DO_SURFACE_EMISSION,    & ! Input flags sources
          DO_PLANE_PARALLEL, DO_2S_LEVELOUT,                             & ! Input flags (including 2p2 7/17/13)
          DO_MVOUT_ONLY, DO_ADDITIONAL_MVOUT, DO_POSTPROCESSING,         & ! Input !@@ New line, 2p3
          NLAYERS, NTOTAL, NBEAMS, N_USER_STREAMS, FOURIER, THREAD,      & ! Input integer control
          FLUX_FACTOR, STREAM_VALUE, X0, USER_STREAMS, PI4,              & ! Input real control
          ALBEDO, BRDF_F_0, BRDF_F, UBRDF_F,                             & ! Input real surface
          THERMAL_BB_INPUT, SURFBB, EMISSIVITY,                          & ! Input real thermal
          DELTAU_VERT, OMEGA_TOTAL, ASYMM_TOTAL, CHAPMAN_FACTORS,        & ! Input real optical
          INITIAL_TRANS, AVERAGE_SECANT, LOCAL_CSZA, LAYER_PIS_CUTOFF,   & ! In/Out
          DELTAU_SLANT, TAUSLANT, SOLAR_BEAM_OPDEP, DO_DIRECTBEAM,       & ! In/Out
          T_DELT_MUBAR, T_DELT_USERM, EMULT_UP, EMULT_DN,                & ! In/Out
          TCOM1, DELTAU_POWER,                                           & ! In/Out
          INTENSITY_F_UP, INTENSITY_F_DN, RADLEVEL_F_UP, RADLEVEL_F_DN,  & ! Output
          FLUXES_TOA, FLUXES_BOA, STATUS_SUB, E_MESSAGE, E_TRACE_1 )       ! Output (modified 2p3, Fluxes)

!  Exception handling

        IF ( STATUS_SUB .NE. 0 ) THEN
          STATUS_EXECUTION = 1
          WRITE(CF,'(I2)')FOURIER
          E_TRACE_2 = 'Error from 2S_FOURIER_MASTER, Fourier # ' &
                        //CF//', Thread # '//wthread
          RETURN
        ENDIF

!  Fourier summation and Convergence examination
!  SS code not included in this version---------------
!  !@@ Alternative Call for Observationsl Geometry case      !@@ 2p1
!  !@@ Call statements expanded to include ALL-LEVEL outputs !@@ 2p2
!  !@@ Convergence skipped for MVOUT_ONLY option             !@@ 2p3

!mick fix 12/17/2013 - fixed logic
        !IF ( DO_MVOUT_ONLY ) then
        IF ( .not. DO_MVOUT_ONLY ) then
          DO IBEAM = 1, NBEAMS
            IF ( DO_USER_OBSGEOMS.and.DO_SOLAR_SOURCES ) THEN
              CALL TWOSTREAM_CONVERGE_OBSGEO &
              ( MAXBEAMS, MAX_USER_STREAMS, MAX_USER_RELAZMS,      & ! Dimensions
                MAX_GEOMETRIES, MAXTHREADS, MAXLAYERS,             & ! Dimensions ! @@ 2p2
                DO_UPWELLING, DO_DNWELLING, DO_2S_LEVELOUT,        & ! Inputs     ! @@ 2p2
                NLAYERS, THREAD, IBEAM, FOURIER, AZMFAC,           & ! Inputs     ! @@ 2p2
                INTENSITY_F_UP,  INTENSITY_F_DN,                   & ! Inputs
                RADLEVEL_F_UP,   RADLEVEL_F_DN,                    & ! Inputs     ! @@ 2p2
                INTENSITY_TOA, INTENSITY_BOA,                      & ! In/Out
                RADLEVEL_UP,   RADLEVEL_DN   )                       ! In/Out     ! @@ 2p2
            ELSE
              CALL TWOSTREAM_CONVERGE &
              ( MAXBEAMS, MAX_USER_STREAMS, MAX_USER_RELAZMS,   & ! Dimensions
                MAX_GEOMETRIES, MAXTHREADS, MAXLAYERS,          & ! Dimensions ! @@ 2p2
                DO_UPWELLING, DO_DNWELLING, DO_2S_LEVELOUT,     & ! Inputs     ! @@ 2p2
                NLAYERS, THREAD, IBEAM, FOURIER,                & ! Inputs     ! @@ 2p2
                N_USER_STREAMS, N_USER_RELAZMS, AZMFAC, UMOFF,  & ! Inputs
                INTENSITY_F_UP,  INTENSITY_F_DN,                & ! Inputs
                RADLEVEL_F_UP,   RADLEVEL_F_DN,                 & ! Inputs     ! @@ 2p2
                INTENSITY_TOA, INTENSITY_BOA,                   & ! In/Out
                RADLEVEL_UP,   RADLEVEL_DN   )                    ! In/Out     ! @@ 2p2
            ENDIF
          END DO
        ENDIF

!  End Fourier loop

      ENDDO

!  Finish

      RETURN
END SUBROUTINE TWOSTREAM_MASTER

!

SUBROUTINE TWOSTREAM_FOURIER_MASTER &
     ( MAXLAYERS, MAXTOTAL, MAXBEAMS, MAX_USER_STREAMS, MAXTHREADS,   & ! Dimensions
       DO_UPWELLING, DO_DNWELLING, DO_BRDF_SURFACE, DO_USER_OBSGEOMS, & ! Input flags control (including Obsgeom)
       DO_SOLAR_SOURCES, DO_THERMAL_EMISSION, DO_SURFACE_EMISSION,    & ! Input flags sources
       DO_PLANE_PARALLEL, DO_2S_LEVELOUT,                             & ! Input flags (including new 7/17/13)
       DO_MVOUT_ONLY, DO_ADDITIONAL_MVOUT, DO_POSTPROCESSING,         & ! Input !@@ New line, 2p3
       NLAYERS, NTOTAL, NBEAMS, N_USER_STREAMS, FOURIER, THREAD,      & ! Input integer control
       FLUX_FACTOR, STREAM_VALUE, X0, USER_STREAMS, PI4,              & ! Input real control
       ALBEDO, BRDF_F_0, BRDF_F, UBRDF_F,                             & ! Input real surface
       THERMAL_BB_INPUT, SURFBB, EMISS,                               & ! Input real thermal
       DELTAU_VERT, OMEGA_TOTAL, ASYMM_TOTAL, CHAPMAN_FACTORS,        & ! Input real optical
       INITIAL_TRANS, AVERAGE_SECANT, LOCAL_CSZA, LAYER_PIS_CUTOFF,   & ! In/Out
       DELTAU_SLANT, TAUSLANT, SOLAR_BEAM_OPDEP, DO_DIRECTBEAM,       & ! In/Out
       T_DELT_MUBAR, T_DELT_USERM, EMULT_UP, EMULT_DN,                & ! In/Out
       TCOM1, DELTAU_POWER,                                           & ! In/Out
       INTENSITY_F_UP, INTENSITY_F_DN, RADLEVEL_F_UP, RADLEVEL_F_DN,  & ! Output
       FLUXES_TOA, FLUXES_BOA, STATUS, MESSAGE, TRACE )                 ! Output (modified 2p3, Fluxes)

      implicit none

!  precision

      INTEGER, PARAMETER :: dp     = KIND( 1.0D0 )

!  input
!  -----

!  Dimensions :
!      MAXTOTAL  = 2 * MAXLAYERS
     
      INTEGER, INTENT(IN)        :: MAXLAYERS, MAXTOTAL
      INTEGER, INTENT(IN)        :: MAXBEAMS, MAX_USER_STREAMS
      INTEGER, INTENT(IN)        :: MAXTHREADS

!  Flags
!     ! @@ Rob Spurr, 17 July 2013, Version 2.2, Levelout flag

      LOGICAL, INTENT(IN)  :: DO_UPWELLING, DO_DNWELLING
      LOGICAL, INTENT(IN)  :: DO_BRDF_SURFACE, DO_PLANE_PARALLEL
      LOGICAL, INTENT(IN)  :: DO_2S_LEVELOUT

!  ** New **. October 2011, Sources control, including thermal

      LOGICAL, INTENT(IN)  :: DO_THERMAL_EMISSION
      LOGICAL, INTENT(IN)  :: DO_SURFACE_EMISSION
      LOGICAL, INTENT(IN)  :: DO_SOLAR_SOURCES

!   !@@ Observational Geometry flag !@@ 2p1

      LOGICAL, INTENT(IN)  :: DO_USER_OBSGEOMS !@@ 2p1

!  !@@ Version 2p3, 11/5/13. Flux output flags, processing flag

      LOGICAL, INTENT(IN)  :: DO_MVOUT_ONLY
      LOGICAL, INTENT(IN)  :: DO_ADDITIONAL_MVOUT
      LOGICAL, INTENT(IN)  :: DO_POSTPROCESSING

!  Numbers

      INTEGER, INTENT(IN)  :: NLAYERS, NTOTAL
      INTEGER, INTENT(IN)  :: NBEAMS, N_USER_STREAMS

!  Input Fourier component number, THREAD number
!        (latter added 11/5/13 for Version 2p3 Flux output)

      INTEGER, INTENT(IN)        :: FOURIER, THREAD

!  Flux factor

      REAL(kind=dp), INTENT(IN)  :: FLUX_FACTOR

!  Stream value

      REAL(kind=dp), INTENT(IN)  :: STREAM_VALUE

!  Geometry

      REAL(kind=dp), INTENT(IN)  :: X0           ( MAXBEAMS )
      REAL(kind=dp), INTENT(IN)  :: USER_STREAMS ( MAX_USER_STREAMS )

!  4pi

      REAL(kind=dp), INTENT(IN)  :: PI4

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

!  ** New **. October 2011. Thermal variables
!  ------------------------------------------

      REAL(kind=dp), INTENT(IN)  :: SURFBB
      REAL(kind=dp), INTENT(IN)  :: THERMAL_BB_INPUT ( 0:MAXLAYERS )
      REAL(kind=dp), INTENT(IN)  :: EMISS

!  Optical properties
!  ------------------

      REAL(kind=dp), INTENT(IN)  :: DELTAU_VERT(MAXLAYERS)
      REAL(kind=dp), INTENT(IN)  :: OMEGA_TOTAL(MAXLAYERS)
      REAL(kind=dp), INTENT(IN)  :: ASYMM_TOTAL(MAXLAYERS)

!  SS flux multiplier, not required in this version
!      REAL(kind=dp) SS_FLUX_MULTIPLIER

!  Output
!  ------

!  User-defined solutions

      REAL(kind=dp), INTENT(OUT) :: INTENSITY_F_UP (MAX_USER_STREAMS,MAXBEAMS)
      REAL(kind=dp), INTENT(OUT) :: INTENSITY_F_DN (MAX_USER_STREAMS,MAXBEAMS)

!  Flux output (already initialized here)
!     ! @@ Rob Spurr, 05 November 2013, Version 2.3 --> Flux Output

     REAL(kind=dp), INTENT(INOUT) :: FLUXES_TOA(MAXBEAMS,2,MAXTHREADS)
     REAL(kind=dp), INTENT(INOUT) :: FLUXES_BOA(MAXBEAMS,2,MAXTHREADS)

!  Fourier-component solutions at ALL levels
!     ! @@ Rob Spurr, 17 July 2013, Version 2.2 --> Optional Output at ALL LEVELS

      REAL(kind=dp), INTENT(OUT) :: RADLEVEL_F_UP (MAX_USER_STREAMS,MAXBEAMS,0:MAXLAYERS)
      REAL(kind=dp), INTENT(OUT) :: RADLEVEL_F_DN (MAX_USER_STREAMS,MAXBEAMS,0:MAXLAYERS)

!  Single scatter solutions, commented out in this streamlined version
!      REAL(kind=dp) INTENSITY_SS_UP (N_GEOMETRIES)
!      REAL(kind=dp) INTENSITY_SS_DN (N_GEOMETRIES)

!  Exception handling

      INTEGER      , INTENT(OUT)  :: STATUS
      CHARACTER*(*), INTENT(OUT)  :: MESSAGE, TRACE

!  Arrays required at the Top level
!  ================================

!  Solar beam pseudo-spherical setup
!  ---------------------------------

!  Chapman factors (from pseudo-spherical geometry)

      REAL(kind=dp), INTENT(IN) :: CHAPMAN_FACTORS ( MAXLAYERS, MAXLAYERS, MAXBEAMS )

!     Last layer to include Particular integral solution
!     Average-secant and initial tramsittance factors for solar beams.
!     Solar beam attenuation

      INTEGER      , INTENT(INOUT) :: LAYER_PIS_CUTOFF ( MAXBEAMS )
      REAL(kind=dp), INTENT(INOUT) :: INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS )
      REAL(kind=dp), INTENT(INOUT) :: AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS )
      REAL(kind=dp), INTENT(INOUT) :: LOCAL_CSZA     ( MAXLAYERS, MAXBEAMS )
      REAL(kind=dp), INTENT(INOUT) :: SOLAR_BEAM_OPDEP ( MAXBEAMS )

!  Derived optical thickness inputs

      REAL(kind=dp), INTENT(INOUT) :: DELTAU_SLANT ( MAXLAYERS, MAXLAYERS, MAXBEAMS )
      REAL(kind=dp), INTENT(INOUT) :: TAUSLANT   ( 0:MAXLAYERS, MAXBEAMS )

!  reflectance flags

      LOGICAL      , INTENT(INOUT) :: DO_DIRECTBEAM ( MAXBEAMS )

!  Transmittance Setups
!  --------------------

!  Transmittance factors for average secant stream
!    Computed in the initial setup stage for Fourier m = 0

      REAL(kind=dp), INTENT(INOUT) :: T_DELT_MUBAR ( MAXLAYERS, MAXBEAMS )

!  Transmittance factors for user-defined stream angles
!    Computed in the initial setup stage for Fourier m = 0

      REAL(kind=dp), INTENT(INOUT) :: T_DELT_USERM ( MAXLAYERS, MAX_USER_STREAMS )

!  Multiplier arrays
!  -----------------

!  Forcing term multipliers (saved for whole atmosphere)

      REAL(kind=dp), INTENT(INOUT) :: EMULT_UP (MAX_USER_STREAMS,MAXLAYERS,MAXBEAMS)
      REAL(kind=dp), INTENT(INOUT) :: EMULT_DN (MAX_USER_STREAMS,MAXLAYERS,MAXBEAMS)

!  Thermal help variables

      REAL(kind=dp), INTENT(INOUT) :: TCOM1 ( MAXLAYERS, 2 )
      REAL(kind=dp), INTENT(INOUT) :: DELTAU_POWER ( MAXLAYERS, 2 )

!  Local Arrays for argument passing
!  =================================

!  Geometry arrays
!  ---------------

!  These just save some Polynomial expansions

      REAL(kind=dp) :: ULP  ( MAX_USER_STREAMS )
      REAL(kind=dp) :: POX  ( MAXBEAMS )
      REAL(kind=dp) :: PX0X ( MAXBEAMS )
      REAL(kind=dp) :: PX11, PXSQ

!  Solar beam pseudo-spherical setup
!  ---------------------------------

!  Atmospheric attenuation

      REAL(kind=dp) :: ATMOS_ATTN ( MAXBEAMS )

!  Direct beam solutions. No USER-term required, MS-mode only

      REAL(kind=dp) :: DIRECT_BEAM ( MAXBEAMS )

!  Transmittance factor

      REAL(kind=dp) :: ITRANS_USERM ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )

!  Multiplier arrays
!  -----------------

!  coefficient functions for user-defined angles

      REAL(kind=dp) :: SIGMA_P(MAXLAYERS,MAX_USER_STREAMS,MAXBEAMS)
      REAL(kind=dp) :: SIGMA_M(MAXLAYERS,MAX_USER_STREAMS,MAXBEAMS)

!  L'Hopital's rule logical variables

      LOGICAL       :: EMULT_HOPRULE (MAXLAYERS,MAX_USER_STREAMS,MAXBEAMS)

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

!  Eigensolutions

      REAL(kind=dp) :: EIGENVALUE(MAXLAYERS)
      REAL(kind=dp) :: EIGENTRANS(MAXLAYERS)

!  Eigenvector solutions

      REAL(kind=dp) :: XPOS(2,MAXLAYERS)
      REAL(kind=dp) :: XNEG(2,MAXLAYERS)

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

!  Beam solution vector, Solutions at layer boundaries

      REAL(kind=dp) :: WVEC(2,MAXLAYERS)
      REAL(kind=dp) :: WUPPER(2,MAXLAYERS)
      REAL(kind=dp) :: WLOWER(2,MAXLAYERS)

!  Auxiliary vectors

      REAL(kind=dp) :: QDIFVEC(MAXLAYERS)
      REAL(kind=dp) :: QSUMVEC(MAXLAYERS)
      REAL(kind=dp) :: QVEC   (MAXLAYERS)

!  Downwelling BOA solution, before reflectance

      REAL(kind=dp) :: H_PARTIC

!  Saved help variables

      REAL(kind=dp) :: W_HELP(0:1)

!  Diffuse-term Particular beam solutions at user-defined angles

      REAL(kind=dp) :: U_WPOS2(MAX_USER_STREAMS,MAXLAYERS)
      REAL(kind=dp) :: U_WNEG2(MAX_USER_STREAMS,MAXLAYERS)

!  Single-scatter Particular beam solutions at user-defined angles
!    ****** NOT REQUIRED for MS-mode only
!      REAL(kind=dp) :: U_WPOS1(MAX_USER_STREAMS,MAXLAYERS)
!      REAL(kind=dp) :: U_WNEG1(MAX_USER_STREAMS,MAXLAYERS)

!  Solution constants of integration, and related quantities

      REAL(kind=dp) :: LCON(MAXLAYERS)
      REAL(kind=dp) :: MCON(MAXLAYERS)
      REAL(kind=dp) :: LCON_XVEC(2,MAXLAYERS)
      REAL(kind=dp) :: MCON_XVEC(2,MAXLAYERS)

!  Thermal solutions (Only required once, do not need Intent(INOUT))
!  -----------------

!  Thermal solution at layer boundaries

      REAL(kind=dp) :: T_WUPPER ( 2, MAXLAYERS )
      REAL(kind=dp) :: T_WLOWER ( 2, MAXLAYERS )

!  Complete layer term solutions

      REAL(kind=dp) :: LAYER_TSUP_UP(MAX_USER_STREAMS,MAXLAYERS)
      REAL(kind=dp) :: LAYER_TSUP_DN(MAX_USER_STREAMS,MAXLAYERS)

!  User solutions

      REAL(kind=dp) :: U_TPOS2 ( MAX_USER_STREAMS, MAXLAYERS, 2 )
      REAL(kind=dp) :: U_TNEG2 ( MAX_USER_STREAMS, MAXLAYERS, 2 )

!  Post-processing variables
!  -------------------------

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

      REAL(kind=dp) :: FLUX_MULTIPLIER
      REAL(kind=dp) :: DELTA_FACTOR
      REAL(kind=dp) :: SURFACE_FACTOR

!  Error tracing

      INTEGER       :: STATUS_SUB

!  ##############
!  initialization
!  ##############

!  Exception handling initialize

      STATUS = 0
      MESSAGE = ' '
      TRACE   = ' '

!  Initialize new output. NOT EFFICIENT - MICK, any suggestions ???
!     ! @@ Rob Spurr, 17 July 2013, Version 2.2 --> Optional output at ALL LEVELS

      RADLEVEL_F_UP = 0.0d0
      RADLEVEL_F_DN = 0.0d0

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
        IF ( FOURIER .EQ. 0 ) DO_INCLUDE_SURFACE = .TRUE.
      ENDIF

!  Direct beam flag (only if above surface flag has been set)

      IF ( DO_INCLUDE_SURFACE ) THEN
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
        SURFACE_FACTOR = 2.0d0
        DELTA_FACTOR   = 1.0d0
      ELSE
        SURFACE_FACTOR = 1.0d0
        DELTA_FACTOR   = 2.0d0
      ENDIF

!  Flux multiplier
!   = 1 / 4.pi with beam sources, 1.0 for thermal

      FLUX_MULTIPLIER   = DELTA_FACTOR

!  ###################################
!  Set up operations (for Fourier = 0)
!  ###################################

      IF ( FOURIER .EQ. 0 ) THEN

!   MISCSETUPS (4 subroutines)  :
!       average-secant formulation,
!       transmittance setup
!       Thermal setup
!       Beam solution multipliers

!  Prepare quasi-spherical attenuation

        CALL TWOSTREAM_QSPREP &
       ( MAXLAYERS, MAXBEAMS,                        & ! Dimensions
         NLAYERS, NBEAMS, DO_PLANE_PARALLEL,         & ! Input
         DELTAU_VERT, CHAPMAN_FACTORS, X0,           & ! Input
         DO_DIRECTBEAM,                              & ! In/Out
         INITIAL_TRANS, AVERAGE_SECANT,              & ! Output
         LOCAL_CSZA, LAYER_PIS_CUTOFF,               & ! Output
         DELTAU_SLANT, TAUSLANT,                     & ! Output
         SOLAR_BEAM_OPDEP )                            ! Output

!  Transmittances and Transmittance factors. !@@ Add flag for Observation Geometry
!    !@@ Add Post-processing flag, 11/5/13

        CALL TWOSTREAM_PREPTRANS &
        ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS,                      & ! Dimensions
          DO_USER_OBSGEOMS, DO_POSTPROCESSING,                        & ! Input flags (2p1,2p3)
          NLAYERS, N_USER_STREAMS, NBEAMS, DELTAU_VERT, USER_STREAMS, & ! Input
          INITIAL_TRANS, AVERAGE_SECANT, LAYER_PIS_CUTOFF,            & ! Input
          T_DELT_MUBAR, T_DELT_USERM, ITRANS_USERM )                    ! Output

!  ** New, October 2011 **. Thermal setup

        IF ( DO_INCLUDE_THERMEMISS ) THEN
          CALL TWOSTREAM_THERMALSETUP &
          ( MAXLAYERS,                                            & ! Dimensions
            NLAYERS, OMEGA_TOTAL, DELTAU_VERT, THERMAL_BB_INPUT,  & ! Input
            DELTAU_POWER, TCOM1 )                                   ! Output
        ENDIF

!   EMULT_MASTER  : Beam source function multipliers.
!      !@@ Add alternative for Observational geometry, 2p1
!      !@@ Avoid altogether if no post-processing

        IF ( DO_SOLAR_SOURCES.and.DO_POSTPROCESSING ) THEN
          IF ( DO_USER_OBSGEOMS ) THEN
            CALL TWOSTREAM_EMULTMASTER_OBSGEO &
            ( MAXLAYERS, MAXBEAMS, MAX_USER_STREAMS,          & ! Dimensions
              DO_UPWELLING, DO_DNWELLING,                     & ! Input
              NLAYERS, NBEAMS, DELTAU_VERT,                   & ! Input
              USER_STREAMS, T_DELT_MUBAR, T_DELT_USERM,       & ! Input
              ITRANS_USERM, AVERAGE_SECANT, LAYER_PIS_CUTOFF, & ! Input
              SIGMA_M, SIGMA_P, EMULT_HOPRULE,                & ! Output
              EMULT_UP, EMULT_DN )                              ! Output
          ELSE
            CALL TWOSTREAM_EMULTMASTER &
           ( MAXLAYERS, MAXBEAMS, MAX_USER_STREAMS,          & ! Dimensions
             DO_UPWELLING, DO_DNWELLING,                     & ! Input
             NLAYERS, NBEAMS, N_USER_STREAMS, DELTAU_VERT,   & ! Input
             USER_STREAMS, T_DELT_MUBAR, T_DELT_USERM,       & ! Input
             ITRANS_USERM, AVERAGE_SECANT, LAYER_PIS_CUTOFF, & ! Input
             SIGMA_M, SIGMA_P, EMULT_HOPRULE,                & ! Output
             EMULT_UP, EMULT_DN )                              ! Output
          ENDIF
        ENDIF

!  End setups operation

      ENDIF

!  debug, 28 dec 12
!      do n = 1, nlayers
!        write(*,*)n,EMULT_DN(1,n,1),EMULT_DN(1,n,2)
!      enddo
!     pause'Emults'

!  Reflected Direct beam attenuation.

      CALL TWOSTREAM_DIRECTBEAM & 
          ( MAXBEAMS,                             & ! Dimensions
            DO_INCLUDE_SURFACE, DO_BRDF_SURFACE,  & ! Input
            NBEAMS, FOURIER, FLUX_FACTOR, X0,     & ! Input
            DELTA_FACTOR, ALBEDO, BRDF_F_0,       & ! Input
            SOLAR_BEAM_OPDEP, DO_DIRECTBEAM,      & ! Input
            ATMOS_ATTN, DIRECT_BEAM )               ! Output

!  Auxiliary Geometry
!  ! @@2p3, 11/5/13/ add Post-processing flag

      CALL TWOSTREAM_AUXGEOM &
      ( MAX_USER_STREAMS, MAXBEAMS, DO_POSTPROCESSING, & ! Dimensions, Flag (2p3
        N_USER_STREAMS, NBEAMS, FOURIER, & ! Input
        X0, USER_STREAMS, STREAM_VALUE,  & ! Input
        PX11, PXSQ, POX, PX0X, ULP )       ! Output

!  #########################################
!   RTE HOMOGENEOUS SOLUTIONS and BVP SETUP
!  #########################################

!  Start layer loop

      DO N = 1, NLAYERS

!  Get Discrete ordinate solutions for this layer

        CALL TWOSTREAM_HOM_SOLUTION &
          ( MAXLAYERS,                               & ! Dimensions
            N, FOURIER, STREAM_VALUE, PXSQ,          & ! Input
            OMEGA_TOTAL, ASYMM_TOTAL, DELTAU_VERT,   & ! Input
            SAB, DAB, EIGENVALUE, EIGENTRANS,        & ! In/Out
            XPOS, XNEG )                               ! In/Out

!  Get Post-processing ("user") solutions for this layer
!   !@@ 2p3. 11/5/13. Post-processing control

        IF ( DO_POSTPROCESSING ) THEN
           CALL TWOSTREAM_HOM_USERSOLUTION &
            ( MAXLAYERS, MAX_USER_STREAMS,                             & ! Dimensions
              N_USER_STREAMS, N, FOURIER, STREAM_VALUE, PX11,          & ! Input
              USER_STREAMS, ULP, XPOS, OMEGA_TOTAL, ASYMM_TOTAL,       & ! Input
              U_XPOS, U_XNEG, &                                          ! In/Out
              U_HELP_P, U_HELP_M )                                       ! Output
        ENDIF

!  end layer loop

      ENDDO

!  Prepare homogeneous solution multipliers
!   !@@ 2p3. 11/5/13. Post-processing control

      IF ( DO_POSTPROCESSING ) THEN
        CALL TWOSTREAM_HMULT_MASTER &
           ( MAXLAYERS, MAX_USER_STREAMS,           & ! Dimensions
             NLAYERS, N_USER_STREAMS, USER_STREAMS, & ! Input
             EIGENVALUE, EIGENTRANS, T_DELT_USERM,  & ! Input
             ZETA_M, ZETA_P, HMULT_1, HMULT_2 )       ! Output
      ENDIF

!  Boundary value problem - MATRIX PREPARATION (Pentadiagonal solution)

      CALL TWOSTREAM_BVP_MATSETUP_PENTADIAG &
         ( MAXLAYERS, MAXTOTAL,                                     & ! Dimensions
           DO_INCLUDE_SURFACE, FOURIER, NLAYERS, NTOTAL,            & ! Input
           DO_BRDF_SURFACE, SURFACE_FACTOR, ALBEDO, BRDF_F,         & ! Input
           XPOS, XNEG, EIGENTRANS, STREAM_VALUE,                    & ! Input
           H_HOMP, H_HOMM, MAT, ELM, SELM,                          & ! Output
           STATUS_SUB, MESSAGE )                                      ! Output

!  Exception handling for Pentadiagonal setup

      IF ( STATUS_SUB .NE. 0 ) THEN
        TRACE  = 'Call BVP_MATSETUP_PENTADIAG in 2S_FOURIER_MASTER'
        STATUS = 1
        RETURN
      ENDIF

!  Thermal solutions
!     1. Find the Particular solution (NOT FOR transmittance only)
!     2. Compute thermal layer source terms. (Upwelling and Downwelling)
!       These will be scaled up by factor 4.pi if solar beams as well

!   !@@ 2p3. 11/5/13. Post-processing control

      IF ( DO_INCLUDE_THERMEMISS ) THEN

        CALL TWOSTREAM_THERMALSOLUTION &
        ( MAXLAYERS, MAX_USER_STREAMS, DO_POSTPROCESSING,       & ! Dimensions
          DO_UPWELLING, DO_DNWELLING, NLAYERS, N_USER_STREAMS,  & ! Input
          STREAM_VALUE, USER_STREAMS, OMEGA_TOTAL, ASYMM_TOTAL, & ! Input
          SAB, DAB, DELTAU_POWER, TCOM1,                        & ! Input
          T_WUPPER, T_WLOWER, U_TPOS2, U_TNEG2 )                  ! Output

        IF ( DO_POSTPROCESSING ) THEN
          CALL TWOSTREAM_THERMALSTERMS &
          ( MAXLAYERS, MAX_USER_STREAMS,                    & ! Dimensions
            DO_UPWELLING, DO_DNWELLING, DO_SOLAR_SOURCES,   & ! Input
            NLAYERS, N_USER_STREAMS, PI4, USER_STREAMS,     & ! Input
            T_DELT_USERM, DELTAU_POWER, U_TPOS2, U_TNEG2,   & ! Input
            LAYER_TSUP_UP, LAYER_TSUP_DN  )                   ! Output
        ENDIF

      ENDIF

!  Skip the thermal-only section if there are solar sources

      IF ( DO_SOLAR_SOURCES ) GO TO 455

!  ####################################################
!  Complete Radiation Field with Thermal-only solutions
!  ####################################################

!  Only one solution, local direct_beam flag NOT set

      IBEAM = 1
      DO_INCLUDE_DIRECTBEAM = .FALSE.

!  set the BVP PI solution at the lower/upper boundaries

      DO N = 1, NLAYERS
         DO I = 1, 2
            WUPPER(I,N) = T_WUPPER(I,N)
            WLOWER(I,N) = T_WLOWER(I,N)
         ENDDO
      ENDDO

!  Solve boundary value problem (Pentadiagonal solution)

      CALL TWOSTREAM_BVP_SOLUTION_PENTADIAG &
      ( MAXLAYERS, MAXBEAMS, MAXTOTAL,                       & ! Dimensions
        DO_INCLUDE_SURFACE, DO_INCLUDE_DIRECTBEAM,           & ! Input
        DO_INCLUDE_SURFEMISS, DO_BRDF_SURFACE,               & ! Input
        FOURIER, IBEAM, NLAYERS, NTOTAL,                     & ! Input
        SURFACE_FACTOR, ALBEDO, BRDF_F, EMISS, SURFBB,       & ! Input
        DIRECT_BEAM, XPOS, XNEG, WUPPER, WLOWER,             & ! Input
        STREAM_VALUE, MAT, ELM, SELM,                        & ! Input
        H_PARTIC, LCON, MCON, LCON_XVEC, MCON_XVEC )           ! Output

!  upwelling, MSMODE only, no Direct Beam inclusion.
!         !@@ 2p1 New OBSGEOM     option 12/21/12
!         !@@ 2p2 New 2S_LEVELOUT option 07/17/13
!         !@@ 2p3. 11/5/13. Post-processing control

      IF ( DO_UPWELLING .and. DO_POSTPROCESSING ) THEN
        CALL TWOSTREAM_UPUSER_INTENSITY &
            ( MAXLAYERS, MAXBEAMS, MAX_USER_STREAMS,                   & ! Dimensions
              DO_INCLUDE_SURFACE, DO_BRDF_SURFACE, DO_USER_OBSGEOMS,   & ! Input !@@ 2p1
              DO_SOLAR_SOURCES, DO_INCLUDE_THERMEMISS, DO_2S_LEVELOUT, & ! Input !@@ 2p2
              FOURIER, IBEAM, NLAYERS, N_USER_STREAMS,                 & ! Input
              SURFACE_FACTOR, ALBEDO, UBRDF_F,                         & ! Input
              FLUX_MULTIPLIER, PI4, T_DELT_USERM, STREAM_VALUE,        & ! Input
              EIGENTRANS, LCON, LCON_XVEC, MCON, MCON_XVEC,            & ! Input
              WLOWER, U_XPOS, U_XNEG, U_WPOS2,                         & ! Input
              HMULT_1, HMULT_2, EMULT_UP, LAYER_TSUP_UP,               & ! Input
              IDOWNSURF, INTENSITY_F_UP, RADLEVEL_F_UP, CUMSOURCE_UP )   ! Output !@@ 2p2
      ENDIF

!  Downwelling, MSMODE only,
!         !@@ 2p1 New OBSGEOM     option 12/21/12
!         !@@ 2p2 New 2S_LEVELOUT option 07/17/13
!         !@@ 2p3. 11/5/13. Post-processing control

      IF ( DO_DNWELLING .and. DO_POSTPROCESSING ) THEN
        CALL TWOSTREAM_DNUSER_INTENSITY &
            ( MAXLAYERS, MAXBEAMS, MAX_USER_STREAMS,                     & ! Dimensions
              DO_INCLUDE_THERMEMISS, DO_SOLAR_SOURCES,                   & ! Dimensions
              DO_USER_OBSGEOMS, DO_2S_LEVELOUT,                          & ! Inputs !@@ 2p1, 2p2
              FOURIER, IBEAM, NLAYERS, N_USER_STREAMS,                   & ! Input
              FLUX_MULTIPLIER, PI4, T_DELT_USERM,                        & ! Input
              LCON, MCON, U_XPOS, U_XNEG, U_WNEG2,                       & ! Input
              HMULT_1, HMULT_2, EMULT_DN, LAYER_TSUP_DN,                 & ! Input
              INTENSITY_F_DN, RADLEVEL_F_DN, CUMSOURCE_DN )                ! Output !@@ 2p2
      ENDIF

!  Flux output. New Subroutine, 11/5/13 Version 2.3

      IF ( DO_INCLUDE_MVOUT ) THEN
        CALL TWOSTREAM_FLUXES &
          ( MAXBEAMS, MAXLAYERS, MAXTHREADS, DO_UPWELLING, DO_DNWELLING, & ! Input  Dimensions, flags
            IBEAM, NLAYERS, THREAD, PI4, STREAM_VALUE, FLUX_MULTIPLIER,  & ! Input Control
            LCON_XVEC, MCON_XVEC, EIGENTRANS, WUPPER, WLOWER,            & ! Input 2-stream solution
            FLUXES_TOA, FLUXES_BOA )                                       ! Output
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

!  stream solution

          CALL TWOSTREAM_BEAM_SOLUTION &
          ( MAXLAYERS, MAXBEAMS,                               & ! Dimensions
            N, FOURIER, IBEAM,                                 & ! Input
            FLUX_FACTOR, LAYER_PIS_CUTOFF, STREAM_VALUE, PX0X, & ! Input
            AVERAGE_SECANT, INITIAL_TRANS, T_DELT_MUBAR,       & ! Input
            OMEGA_TOTAL, ASYMM_TOTAL, SAB, DAB, EIGENVALUE,    & ! Input
            QSUMVEC, QDIFVEC, QVEC, WVEC, WUPPER, WLOWER )       ! In/Out

!  user solutions. !@@ 2p1, Additional option for Observation Geometry, 12/21/12
!         !@@ 2p3. 11/5/13. Post-processing control

          IF (  DO_POSTPROCESSING ) THEN
            CALL TWOSTREAM_BEAM_USERSOLUTION &
            ( MAXLAYERS, MAXBEAMS, MAX_USER_STREAMS,             & ! Dimensions
              DO_UPWELLING, DO_DNWELLING, DO_USER_OBSGEOMS,      & ! Input !@@
              N_USER_STREAMS, N, FOURIER, IBEAM,                 & ! Input
              FLUX_FACTOR, LAYER_PIS_CUTOFF, STREAM_VALUE, PX11, & ! Input
              OMEGA_TOTAL, ASYMM_TOTAL, USER_STREAMS, ULP, WVEC, & ! Input
              U_WPOS2, U_WNEG2, &                                  ! In/Out
              W_HELP )                                             ! Output
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

!  Solve boundary value problem (Pentadiagonal solution)
!  ----------------------------

        CALL TWOSTREAM_BVP_SOLUTION_PENTADIAG &
      ( MAXLAYERS, MAXBEAMS, MAXTOTAL,                       & ! Dimensions
        DO_INCLUDE_SURFACE, DO_DIRECTBEAM(IBEAM),            & ! Input
        DO_INCLUDE_SURFEMISS, DO_BRDF_SURFACE,               & ! Input
        FOURIER, IBEAM, NLAYERS, NTOTAL,                     & ! Input
        SURFACE_FACTOR, ALBEDO, BRDF_F, EMISS, SURFBB,       & ! Input
        DIRECT_BEAM, XPOS, XNEG, WUPPER, WLOWER,             & ! Input
        STREAM_VALUE, MAT, ELM, SELM,                        & ! Input
        H_PARTIC, LCON, MCON, LCON_XVEC, MCON_XVEC )           ! Output

! ##################################
!   Radiance Field Post Processing
! ##################################

!  upwelling, MSMODE only, no Direct Beam inclusion.
!         !@@ 2p1 New OBSGEOM     option 12/21/12
!         !@@ 2p2 New 2S_LEVELOUT option 07/17/13
!         !@@ 2p3. 11/5/13. Post-processing control

        IF ( DO_UPWELLING .and. DO_POSTPROCESSING ) THEN
          CALL TWOSTREAM_UPUSER_INTENSITY &
            ( MAXLAYERS, MAXBEAMS, MAX_USER_STREAMS,                   & ! Dimensions
              DO_INCLUDE_SURFACE, DO_BRDF_SURFACE, DO_USER_OBSGEOMS,   & ! Input !@@ 2p1
              DO_SOLAR_SOURCES, DO_INCLUDE_THERMEMISS, DO_2S_LEVELOUT, & ! Input !@@ 2p2
              FOURIER, IBEAM, NLAYERS, N_USER_STREAMS,                 & ! Input
              SURFACE_FACTOR, ALBEDO, UBRDF_F,                         & ! Input
              FLUX_MULTIPLIER, PI4, T_DELT_USERM, STREAM_VALUE,        & ! Input
              EIGENTRANS, LCON, LCON_XVEC, MCON, MCON_XVEC,            & ! Input
              WLOWER, U_XPOS, U_XNEG, U_WPOS2,                         & ! Input
              HMULT_1, HMULT_2, EMULT_UP, LAYER_TSUP_UP,               & ! Input
              IDOWNSURF, INTENSITY_F_UP, RADLEVEL_F_UP, CUMSOURCE_UP )   ! Output !@@ 2p2
        ENDIF

!  Downwelling, MSMODE only,
!         !@@ 2p1 New OBSGEOM     option 12/21/12
!         !@@ 2p2 New 2S_LEVELOUT option 07/17/13
!         !@@ 2p3. 11/5/13. Post-processing control

        IF ( DO_DNWELLING .and. DO_POSTPROCESSING ) THEN
          CALL TWOSTREAM_DNUSER_INTENSITY &
            ( MAXLAYERS, MAXBEAMS, MAX_USER_STREAMS,                     & ! Dimensions
              DO_INCLUDE_THERMEMISS, DO_SOLAR_SOURCES,                   & ! Dimensions
              DO_USER_OBSGEOMS, DO_2S_LEVELOUT,                          & ! Inputs !@@ 2p1, 2p2
              FOURIER, IBEAM, NLAYERS, N_USER_STREAMS,                   & ! Input
              FLUX_MULTIPLIER, PI4, T_DELT_USERM,                        & ! Input
              LCON, MCON, U_XPOS, U_XNEG, U_WNEG2,                       & ! Input
              HMULT_1, HMULT_2, EMULT_DN, LAYER_TSUP_DN,                 & ! Input
              INTENSITY_F_DN, RADLEVEL_F_DN, CUMSOURCE_DN )                ! Output !@@ 2p2
        ENDIF

!  Flux output. New Subroutine, 11/5/13 Version 2.3

        IF ( DO_INCLUDE_MVOUT ) THEN
          CALL TWOSTREAM_FLUXES &
            ( MAXBEAMS, MAXLAYERS, MAXTHREADS, DO_UPWELLING, DO_DNWELLING, & ! Input  Dimensions, flags
              IBEAM, NLAYERS, THREAD, PI4, STREAM_VALUE, FLUX_MULTIPLIER,  & ! Input Control
              LCON_XVEC, MCON_XVEC, EIGENTRANS, WUPPER, WLOWER,            & ! Input 2-stream solution
              FLUXES_TOA, FLUXES_BOA )                                       ! Output
        ENDIF

!  End loop over beam solutions

      END DO

!  ######
!  finish
!  ######

      RETURN
END SUBROUTINE TWOSTREAM_FOURIER_MASTER

end module twostream_master_m
