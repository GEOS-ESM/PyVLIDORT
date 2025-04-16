module 2STREAM_Mod

      USE twostream_pars

      implicit NONE

      PUBLIC TWOSTREAM_Init

      TYPE TWOSTREAM_IO
!  Directional Flags

      LOGICAL       :: DO_UPWELLING
      LOGICAL       :: DO_DNWELLING

!  Plane parallel and deltam-2stream scaling flags

      LOGICAL       :: DO_PLANE_PARALLEL
      LOGICAL       :: DO_D2S_SCALING

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

!  BVP control --- New 6/25/14, Version 2.3 and higher
!   * PentaDiagonal Inverse flag (BVP solved from bottom to top). Only for BVPIndex = 1
!   * BVP Index : 0 = LAPACK, 1 = Penta # 1 (original), 2 = Penta # 2 (new, 2012 Kanal paper)
!   * BVP Scale Factor. Debug only. Set this to 1.0 on input

      LOGICAL        :: DO_PENTADIAG_INVERSE
      INTEGER        :: BVPINDEX
      REAL(kind=dp)  :: BVPSCALEFACTOR

!  Numbers

      INTEGER       :: NLAYERS, NTOTAL

!  Stream value

      REAL(kind=dp) :: STREAM_VALUE

!  Thermal inputs

      REAL(kind=dp) :: SURFBB
      REAL(kind=dp) :: THERMAL_BB_INPUT  ( 0:MAXLAYERS )

!  Emissivity

      REAL(kind=dp) :: EMISSIVITY

!  Version 2p3. 1/23/14. Introduce SLEAVE stuff
!    Do not require any first-order inputs (exact or Fourier)
!    Isotropic Surface leaving term (if flag set)
!  Fourier components of Surface-leaving terms:

      REAL(kind=dp) ::  SLTERM_ISOTROPIC ( MAXBEAMS )
      REAL(kind=dp) ::  SLTERM_F_0 ( 0:1, MAXBEAMS )


!  Flux factor

      REAL(kind=dp) :: FLUX_FACTOR

!  Earth radius

      REAL(kind=dp) :: EARTH_RADIUS

!  Exception handling

!    1. Up to 100 Check Messages and actions

      INTEGER       :: STATUS_INPUTCHECK
      INTEGER       :: C_NMESSAGES
      CHARACTER (LEN=100) :: C_MESSAGES ( 0:MAXMESSAGES )
      CHARACTER (LEN=100) :: C_ACTIONS  ( 0:MAXMESSAGES )

!    2. Execution message and 2 Traces

      INTEGER       :: STATUS_EXECUTION
      CHARACTER (LEN=100) :: E_MESSAGE, E_TRACE_1, E_TRACE_2

      END TYPE TWOSTREAM_IO

      !  2STREAM input structure
      TYPE TWOSTREAM
        logical     :: initialized = .false.
        integer     :: NBEAMS = 1          ! Number of solar zenith angles
        integer     :: N_USER_ANGLES = 1   ! Number of Viewing zenith angles
        integer     :: N_USER_RELAZMS = 1  ! Number of relative azimuth angles
        integer     :: N_USER_LEVELS  = 1  ! Number of user-defined vertical output levels
        integer     :: N_USER_OBSGEOMS = 1 ! Number of observation geometry triplets

        TYPE(TWOSTREAM_IO)      :: TSIO
      END TYPE TWOSTREAM

!.............................................................................

      subroutine TWOSTREAM_Init (self, km, rc)

      USE TWOSTREAM_PARS

      type(TWOSTREAM),              intent(inout) :: self

      integer,                      intent(in)    :: km   ! number of atmospheric layers
      integer,                      intent(out)   :: rc     ! error code

      self%initialized = .true.


!    Customize those parameters that do not vary from pixel to pixel
!     ---------------------------------------------------------------


!                      User-defined output control
!                      ---------------------------

      self%TSIO%DO_UPWELLING     = .true.     ! Upwelling output?
      self%TSIO%DO_DNWELLING     = .false.    ! Downwelling output?
      self%TSIO%DO_USER_OBSGEOMS = .true.     ! Do Observation Geometry?
      self%TSIO%DO_2S_LEVELOUT   = .false.    ! Output at all layer boundaries

      self%TSIO%DO_ADDITIONAL_MVOUT = .true.  ! output actinic and regular fluxes 
                                    ! (upwelling at the TOA and downwelling at the BOA) additionally?
      self%TSIO%DO_MVOUT_ONLY       = .false. ! output actinic and regular fluxes only?

!                            Solar Sources
!                            -------------

      self%TSIO%DO_SOLAR_SOURCES    = .true.     ! Include solar sources?
      self%TSIO%DO_PLANE_PARALLEL   = .false.    ! Plane-parallel treatment of direct beam?

!                          Thermal controls
!                          ----------------

      self%TSIO%DO_THERMAL_EMISSION  = .false.  ! Do thermal emission?
      self%TSIO%DO_SURFACE_EMISSION  =  .false. ! Do Surface emission?

!                         BRDF controls
!                         --------------
      self%TSIO%DO_BRDF_SURFACE    = .false.    ! required to be set here, but this is overwritten in 2STREAM_SurfaceMod.F90
      self%TSIO%DO_SURFACE_LEAVING = .false.    ! required to be set here, but this is overwritten in 2STREAM_SurfaceMod.F90
      self%TSIO%DO_SL_ISOTROPIC    = .false.    ! required to be set here, but this is overwritten in 2STREAM_SurfaceMod.F90

!                         Computational Input Controls
!                         ----------------------------
      self%TSIO%STREAM_VALUE  = 0.5d0  ! Value of stream angle used for the 2S computations. 
                             ! For scattering applications, should be set to 0.5; 
                             ! for thermal applications, should be set to sqrt(1/3)

!                         Performance Control
!                         -------------------
      self%TSIO%DO_D2S_SCALING       = .true.  ! Include Delta-M scaling?
      self%TSIO%DO_PENTADIAG_INVERSE = .false. ! Boundary value problem method test variable.

!  Atmosphere
      self%TSIO%EARTH_RADIUS = 6371.0d0         ! Earth radius (km)
      self%TSIO%NLAYERS      = km               ! Number of atmospheric layers
      self%TSIO%NTOTAL       = 2 * NLAYERS

! Flux factor
      self%TSIO%FLUX_FACTOR  = 1.00d0           ! Solar flux constant; =1 if no solar sources.

! Taylor Control
      self%TSIO%TAYLOR_ORDER = 3                ! Number of small-number terms in Taylor series expansions
    
! BVP value problem method test variables
      self%TSIO%BVPINDEX       = 1              ! Boundary value problem method test variable.
      self%TSIO%BVPSCALEFACTOR = 1.0d0          ! Boundary value problem method test variable.

! Initialize Thermal inputs
      self%TSIO%THERMAL_BB_INPUT = 0.0d0
      self%TSIO%EMISSIVITY       = 0.0d0
      self%TSIO%SURFBB           = 0.0d0

! Initialize SLEAVE inputs
      self%TSIO%SLTERM_ISOTROPIC = 0.0d0
      self%TSIO%SLTERM_F_0       = 0.0d0

      end subroutine TWOSTREAM_Init


      END MODULE twostream_init
