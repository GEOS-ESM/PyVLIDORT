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
! ###########################################################

! #############################################################
! #                                                           #
! #   This Version of LIDORT-2STREAM comes with a GNU-style   #
! #   license. Please read the license carefully.             #
! #                                                           #
! #############################################################

! ###########################################################
! #                                                         #
! #   Contains the following Master subroutines             #
! #                                                         #
! #          TWOSTREAM_UPUSER_INTENSITY   (master)          #
! #          TWOSTREAM_DNUSER_INTENSITY   (master)          #
! #          TWOSTREAM_FLUXES . 11/5/13 Version 2p3         #
! #          TWOSTREAM_CONVERGE (master)                    #
! #          TWOSTREAM_CONVERGE_OBSGEO (master) !@@ 2p1     #
! #                                                         #
! ###########################################################

module twostream_intensity_m

PUBLIC

contains

SUBROUTINE TWOSTREAM_UPUSER_INTENSITY &
  ( MAXLAYERS, MAXBEAMS, MAX_USER_STREAMS,                       & ! Dimensions
    DO_INCLUDE_SURFACE, DO_BRDF_SURFACE, DO_USER_OBSGEOMS,       & ! inputs !@@ 2p1
    DO_SOLAR_SOURCES, DO_INCLUDE_THERMEMISS, DO_2S_LEVELOUT,     & ! inputs !@@ 2p2
    FOURIER_COMPONENT, IPARTIC, NLAYERS,                         & ! inputs
    N_USER_STREAMS, SURFACE_FACTOR, ALBEDO, UBRDF_F,             & ! inputs
    FLUX_MULTIPLIER, PI4, T_DELT_USERM, STREAM_VALUE,            & ! inputs
    T_DELT_EIGEN, LCON, LCON_XVEC, MCON, MCON_XVEC,              & ! inputs
    WLOWER, U_XPOS, U_XNEG, U_WPOS2,                             & ! inputs
    HMULT_1, HMULT_2, EMULT_UP, LAYER_TSUP_UP,                   & ! inputs
    IDOWNSURF, INTENSITY_F_UP, RADLEVEL_F_UP, CUMSOURCE_UP )       ! Output !@@ 2p2

      implicit none

!  precision

      INTEGER, PARAMETER :: dp     = KIND( 1.0D0 )

!  Subroutine input arguments
!  --------------------------

!  Dimensions

      INTEGER, INTENT(IN)        :: MAXLAYERS, MAXBEAMS, MAX_USER_STREAMS

!  local surface control flags

      LOGICAL, INTENT(IN)        :: DO_INCLUDE_SURFACE
      LOGICAL, INTENT(IN)        :: DO_BRDF_SURFACE

!   !@@ Observational Geometry flag !@@ 2p1

      LOGICAL, INTENT(IN)        :: DO_USER_OBSGEOMS !@@ 2p1

!  Local source flags

      LOGICAL, INTENT(IN)        :: DO_SOLAR_SOURCES
      LOGICAL, INTENT(IN)        :: DO_INCLUDE_THERMEMISS

!     ! @@ Rob Spurr, 17 July 2013, Version 2.2, Levelout flag

      LOGICAL, INTENT(IN)        :: DO_2S_LEVELOUT

!  Fourier component, beam index

      INTEGER, INTENT(IN)        :: FOURIER_COMPONENT
      INTEGER, INTENT(IN)        :: IPARTIC

!  Numbers

      INTEGER, INTENT(IN)        :: NLAYERS, N_USER_STREAMS

!  Surface stuff

      REAL(kind=dp), INTENT(IN)  :: SURFACE_FACTOR, ALBEDO
      REAL(kind=dp), INTENT(IN)  :: UBRDF_F ( 0:1, MAX_USER_STREAMS )

!  multiplier, 4pi

      REAL(kind=dp), INTENT(IN)  :: FLUX_MULTIPLIER, PI4

!  Transmittance factors for user-defined stream angles

      REAL(kind=dp), INTENT(IN)  :: T_DELT_USERM ( MAXLAYERS, MAX_USER_STREAMS )

!  Stream value

      REAL(kind=dp), INTENT(IN)  :: STREAM_VALUE

!  No USER_DIRECT_BEAM (MSMODE only ===> No Direct BOA source term)
!      DOUBLE PRECISION USER_DIRECT_BEAM ( MAX_USER_STREAMS, MAXBEAMS )

!  transmittance factors for +/- eigenvalues

      REAL(kind=dp), INTENT(IN)  :: T_DELT_EIGEN(MAXLAYERS)

!  Solution constants of integration

      REAL(kind=dp), INTENT(IN)  :: LCON(MAXLAYERS)
      REAL(kind=dp), INTENT(IN)  :: MCON(MAXLAYERS)

!  Solution constants of integration multiplied by homogeneous solutions

      REAL(kind=dp), INTENT(IN)  :: LCON_XVEC(2,MAXLAYERS)
      REAL(kind=dp), INTENT(IN)  :: MCON_XVEC(2,MAXLAYERS)

!  General beam solutions at the Upper/Lower boundary

      REAL(kind=dp), INTENT(IN)  :: WLOWER(2,MAXLAYERS)

!  Eigenvectors defined at user-defined stream angles
!     EP for the positive KEIGEN values, EM for -ve KEIGEN

      REAL(kind=dp), INTENT(IN)  :: U_XPOS(MAX_USER_STREAMS,MAXLAYERS)
      REAL(kind=dp), INTENT(IN)  :: U_XNEG(MAX_USER_STREAMS,MAXLAYERS)

!  Diffuse-term Particular beam solution at user-defined angles

      REAL(kind=dp), INTENT(IN)  :: U_WPOS2(MAX_USER_STREAMS,MAXLAYERS)

!  Single-scatter Particular beam solution at user-defined angles
!    @@@ NOT REQUIRED For MS-mode only
!     REAL(kind=dp), INTENT(IN)  :: U_WPOS1(MAX_USER_STREAMS,MAXLAYERS)

!  solution multipliers

      REAL(kind=dp), INTENT(IN)  :: HMULT_1 (MAX_USER_STREAMS,MAXLAYERS)
      REAL(kind=dp), INTENT(IN)  :: HMULT_2 (MAX_USER_STREAMS,MAXLAYERS)
      REAL(kind=dp), INTENT(IN)  :: EMULT_UP(MAX_USER_STREAMS,MAXLAYERS,MAXBEAMS)

!  Thermal layer source term

      REAL(kind=dp), INTENT(IN)  :: LAYER_TSUP_UP(MAX_USER_STREAMS,MAXLAYERS)

!  Outputs
!  -------

!  Reflectance integrand  a(j).x(j).I(-j)

      REAL(kind=dp), INTENT(OUT) :: IDOWNSURF

!  User-defined solutions
! mick fix 11/7/2012 - change to "inout"
      REAL(kind=dp), INTENT(INOUT) :: INTENSITY_F_UP(MAX_USER_STREAMS,MAXBEAMS)

!  Fourier-component solutions at ALL levels
!     ! @@ Rob Spurr, 17 July 2013, Version 2.2 --> Optional Output at ALL LEVELS

      REAL(kind=dp), INTENT(INOUT) :: RADLEVEL_F_UP (MAX_USER_STREAMS,MAXBEAMS,0:MAXLAYERS)

!  Cumulative source terms

      REAL(kind=dp), INTENT(OUT) :: CUMSOURCE_UP(MAX_USER_STREAMS,0:MAXLAYERS)

!  local variables
!  ---------------

!  help variables

      LOGICAL       :: LUOGSS
      INTEGER       :: UM, N, NC, M, LUM, N1
      REAL(kind=dp) :: LAYERSOURCE ( MAX_USER_STREAMS )
      REAL(kind=dp) :: BOA_SOURCE  ( MAX_USER_STREAMS )
      REAL(kind=dp) :: SHOM, SFOR2, PAR, HOM, KMULT, TM

!  Local user index !@@ Observation Geometry choice, 12/21/12
!  Local user flag  !@@ Observation Geometry choice, 12/21/12

      LUM    = 1
      LUOGSS = DO_USER_OBSGEOMS .and. DO_SOLAR_SOURCES

!  Dummy (for debug)

      M = FOURIER_COMPONENT

!  Zero all Fourier components - New rule, better for safety
!    Only did this for components close to zenith (formerly)
!            !@@ 2p1, Observation Geometry choice, 12/21/12
!            !@@ 2p2, Zero the new "All-level" output (Already zeroed in Fourier Master)

      IF ( LUOGSS ) THEN
        INTENSITY_F_UP(LUM,IPARTIC) = 0.0d0
     ELSE
        DO UM = 1, N_USER_STREAMS
          INTENSITY_F_UP(UM,IPARTIC) = 0.0d0
        ENDDO
      ENDIF

!  BOA source terms
!  ----------------

!  initialise boa source terms
!    MSMODE only ===> No Direct BOA source term
!            !@@ Observation Geometry choice, 12/21/12

      IF ( LUOGSS ) THEN
         BOA_SOURCE(IPARTIC)    = 0.0d0
      ELSE
        DO UM = 1, N_USER_STREAMS
          BOA_SOURCE(UM)        = 0.0d0
        ENDDO
      ENDIF

!  Full solution: Downward intensity at computational angles (beam/homog)
!     --> Develop reflectance integrand  a(j).x(j).I(-j)

      N = NLAYERS
      PAR = WLOWER(1,N)
      HOM = LCON_XVEC(1,N)*T_DELT_EIGEN(N) + MCON_XVEC(1,N)
      IDOWNSURF = ( PAR + HOM ) * STREAM_VALUE

!  BOA source terms
!  ----------------

      IF ( DO_INCLUDE_SURFACE )THEN

!  reflected multiple scatter intensity at user defined-angles
!    BRDF code added 4 May 2009
!            !@@ Observation Geometry choice, 12/21/12

        IF ( DO_BRDF_SURFACE  ) THEN
          KMULT = SURFACE_FACTOR * IDOWNSURF
          IF ( LUOGSS ) THEN
              BOA_SOURCE(IPARTIC) = KMULT * UBRDF_F(M,IPARTIC)
          ELSE
            DO UM = 1, N_USER_STREAMS
              BOA_SOURCE(UM) = KMULT * UBRDF_F(M,UM)
            ENDDO
          ENDIF
        ELSE
          KMULT  = SURFACE_FACTOR * ALBEDO * IDOWNSURF
          IF ( LUOGSS ) THEN
            BOA_SOURCE(IPARTIC) = KMULT
          ELSE
            DO UM = 1, N_USER_STREAMS
              BOA_SOURCE(UM) = KMULT
            ENDDO
          ENDIF
        ENDIF

!    MSMODE only ===> No Direct BOA source term
!        IF ( DO_INCLUDE_DIRECTBEAM ) THEN
!          DO UM = 1, N_USER_STREAMS
!            DIRECT_BOA_SOURCE(UM) = USER_DIRECT_BEAM(UM,IPARTIC)
!          ENDDO
!        ENDIF

      ENDIF

!  Add surface emission term if flagged
!    ********  Direct surface emission not included

!  Initialize post-processing recursion
!  ====================================

!  Set the cumulative source term equal to BOA values
!    MSMODE only ===> No Direct BOA source term
!            !@@ 2p1, Observation Geometry choice, 12/21/12
!            !@@ 2p2, Set All-level output at surface

      NC = 0
      IF ( LUOGSS ) THEN
        CUMSOURCE_UP(IPARTIC,NC) = BOA_SOURCE(IPARTIC)
        if ( DO_2S_LEVELOUT ) RADLEVEL_F_UP(LUM,IPARTIC,NLAYERS) = FLUX_MULTIPLIER * CUMSOURCE_UP(IPARTIC,NC)
      ELSE
        DO UM = 1, N_USER_STREAMS
          CUMSOURCE_UP(UM,NC) = BOA_SOURCE(UM)
          if ( DO_2S_LEVELOUT ) RADLEVEL_F_UP(UM,IPARTIC,NLAYERS) = FLUX_MULTIPLIER * CUMSOURCE_UP(UM,NC)
        ENDDO
      ENDIF

!  Recursion Loop in Source function integration
!  =============================================

      DO N = NLAYERS, 1, -1
        NC = NLAYERS + 1 - N ; N1 = N - 1

!  Homogeneous solutions. !@@ 2p1, Observation Geometry choice, 12/21/12

        IF ( LUOGSS ) THEN
          SHOM = LCON(N) * U_XPOS(IPARTIC,N) * HMULT_2(IPARTIC,N) + &
                 MCON(N) * U_XNEG(IPARTIC,N) * HMULT_1(IPARTIC,N)
          LAYERSOURCE(IPARTIC) = SHOM
        ELSE
          DO UM = 1, N_USER_STREAMS
            SHOM = LCON(N) * U_XPOS(UM,N) * HMULT_2(UM,N) + &
                   MCON(N) * U_XNEG(UM,N) * HMULT_1(UM,N)
            LAYERSOURCE(UM) = SHOM
          ENDDO
        ENDIF

!        if (ipartic.eq.2)write(*,*)n,nc,LAYERSOURCE(IPARTIC)
!        if (n.eq.1.and.ipartic.eq.2)pause

!  Add thermal emission term (direct and diffuse)
!     Modulus 4.pi if solar sources are included (taken care of earlier)
!       !@@ 2p1, Observation Geometry choice, 12/21/12

        IF ( DO_INCLUDE_THERMEMISS ) THEN
          TM = 1.0_dp ; IF ( DO_SOLAR_SOURCES ) TM = 1.0_dp/PI4
          IF ( LUOGSS ) THEN
            LAYERSOURCE(IPARTIC) = LAYERSOURCE(IPARTIC) + LAYER_TSUP_UP(IPARTIC,N)*TM
          ELSE
            DO UM = 1 , N_USER_STREAMS
              LAYERSOURCE(UM) = LAYERSOURCE(UM) + LAYER_TSUP_UP(UM,N)*TM
            ENDDO
          ENDIF
        ENDIF

!  Add solar source term
!       !@@ 2p1, Observation Geometry choice, 12/21/12

        IF ( DO_SOLAR_SOURCES ) THEN
          IF ( DO_USER_OBSGEOMS ) THEN
            SFOR2 =  EMULT_UP(LUM,N,IPARTIC) * U_WPOS2(IPARTIC,N)
            LAYERSOURCE(IPARTIC) = LAYERSOURCE(IPARTIC) + SFOR2
          ELSE
            DO UM = 1, N_USER_STREAMS
              SFOR2 =  EMULT_UP(UM,N,IPARTIC) * U_WPOS2(UM,N)
              LAYERSOURCE(UM) = LAYERSOURCE(UM) + SFOR2
            ENDDO
          ENDIF
        ENDIF

!   This is not required-----------------------No SS correction
!        IF ( .NOT.DO_MSMODE_LIDORT ) THEN
!          DO UM = 1, N_USER_STREAMS
!            SFOR1 = U_WPOS1(UM,N) * EMULT_UP(UM,N,IPARTIC)
!            LAYERSOURCE(UM) = LAYERSOURCE(UM) + SFOR1
!          ENDDO
!        ENDIF

!  Upgrade source function. 
!        !@@ 2p1, Observation Geometry choice, 12/21/12
!        !@@ 2p2, Set All-level outputs, 7/17/13

        IF ( LUOGSS ) THEN
          CUMSOURCE_UP(IPARTIC,NC) = LAYERSOURCE(IPARTIC) + &
                T_DELT_USERM(N,IPARTIC)*CUMSOURCE_UP(IPARTIC,NC-1)
          if (DO_2S_LEVELOUT)RADLEVEL_F_UP(LUM,IPARTIC,N1) = FLUX_MULTIPLIER * CUMSOURCE_UP(IPARTIC,NC)
        ELSE
          DO UM = 1, N_USER_STREAMS
            CUMSOURCE_UP(UM,NC) = LAYERSOURCE(UM) + &
                T_DELT_USERM(N,UM)*CUMSOURCE_UP(UM,NC-1)
            if (DO_2S_LEVELOUT)RADLEVEL_F_UP(UM,IPARTIC,N1) = FLUX_MULTIPLIER * CUMSOURCE_UP(UM,NC)
          ENDDO
        ENDIF

      ENDDO

!  User-defined stream output, just set to the cumulative source term
!    !@@ 2p1, Observation Geometry choice, 12/21/12

      IF ( LUOGSS ) THEN
        INTENSITY_F_UP(LUM,IPARTIC) = FLUX_MULTIPLIER * CUMSOURCE_UP(IPARTIC,NC)
      ELSE
        DO UM = 1, N_USER_STREAMS
          INTENSITY_F_UP(UM,IPARTIC) = FLUX_MULTIPLIER * CUMSOURCE_UP(UM,NC)
        ENDDO
      ENDIF

!  debug 28 dec 12
!      write(*,*)FOURIER_COMPONENT,INTENSITY_F_UP(1,IPARTIC)
!      if ( ipartic.eq.2.and.FOURIER_COMPONENT.gt.0)pause

!  Finish

      RETURN
END SUBROUTINE TWOSTREAM_UPUSER_INTENSITY

!

SUBROUTINE TWOSTREAM_DNUSER_INTENSITY &
   ( MAXLAYERS, MAXBEAMS, MAX_USER_STREAMS,                     & ! Dimensions
     DO_INCLUDE_THERMEMISS, DO_SOLAR_SOURCES,                   & ! Dimensions
     DO_USER_OBSGEOMS, DO_2S_LEVELOUT,                          & ! Inputs !@@ 2p1, 2p2
     FOURIER_COMPONENT, IPARTIC, NLAYERS,                       & ! Inputs
     N_USER_STREAMS, FLUX_MULTIPLIER, PI4,                      & ! Inputs
     T_DELT_USERM, LCON, MCON, U_XPOS, U_XNEG, U_WNEG2,         & ! Inputs
     HMULT_1, HMULT_2, EMULT_DN, LAYER_TSUP_DN,                 & ! Inputs
     INTENSITY_F_DN, RADLEVEL_F_DN, CUMSOURCE_DN )                ! Output !@@ 2p2

      implicit none

!  precision

      INTEGER, PARAMETER :: dp     = KIND( 1.0D0 )

!  Subroutine input arguments
!  --------------------------

!  Dimensions

      INTEGER, INTENT(IN)        :: MAXLAYERS, MAXBEAMS, MAX_USER_STREAMS

!  Local source flags

      LOGICAL, INTENT(IN)        :: DO_SOLAR_SOURCES
      LOGICAL, INTENT(IN)        :: DO_INCLUDE_THERMEMISS

!   !@@ 2p1, Observational Geometry flag

      LOGICAL, INTENT(IN)        :: DO_USER_OBSGEOMS !@@ 2p1

!     ! @@ Rob Spurr, 17 July 2013, Version 2.2, Levelout flag

      LOGICAL, INTENT(IN)        :: DO_2S_LEVELOUT

!  Fourier component, beam index

      INTEGER, INTENT(IN)        :: FOURIER_COMPONENT
      INTEGER, INTENT(IN)        :: IPARTIC

!  Numbers

      INTEGER, INTENT(IN)        :: NLAYERS, N_USER_STREAMS

!  multiplier, 4pi

      REAL(kind=dp), INTENT(IN)  :: FLUX_MULTIPLIER, PI4

!  Transmittance factors for user-defined stream angles

      REAL(kind=dp), INTENT(IN)  :: T_DELT_USERM ( MAXLAYERS, MAX_USER_STREAMS )

!  Solution constants of integration

      REAL(kind=dp), INTENT(IN)  :: LCON(MAXLAYERS)
      REAL(kind=dp), INTENT(IN)  :: MCON(MAXLAYERS)

!  Eigenvectors defined at user-defined stream angles

      REAL(kind=dp), INTENT(IN)  :: U_XPOS(MAX_USER_STREAMS,MAXLAYERS)
      REAL(kind=dp), INTENT(IN)  :: U_XNEG(MAX_USER_STREAMS,MAXLAYERS)

!  Diffuse-term Particular beam solution at user-defined angles

      REAL(kind=dp), INTENT(IN)  :: U_WNEG2(MAX_USER_STREAMS,MAXLAYERS)

!  Single-scatter Particular beam solution at user-defined angles
!    @@@ NOT REQUIRED For MS-mode only
!      REAL(kind=dp), INTENT(IN)  :: U_WNEG1(MAX_USER_STREAMS,MAXLAYERS)

!  solution multipliers 

      REAL(kind=dp), INTENT(IN)  :: HMULT_1(MAX_USER_STREAMS,MAXLAYERS)
      REAL(kind=dp), INTENT(IN)  :: HMULT_2(MAX_USER_STREAMS,MAXLAYERS)
      REAL(kind=dp), INTENT(IN)  :: EMULT_DN(MAX_USER_STREAMS,MAXLAYERS,MAXBEAMS)

!  Thermal layer source term

      REAL(kind=dp), INTENT(IN)  :: LAYER_TSUP_DN(MAX_USER_STREAMS,MAXLAYERS)

!  Outputs
!  -------

!  User-defined solutions

      REAL(kind=dp), INTENT(INOUT) :: INTENSITY_F_DN(MAX_USER_STREAMS,MAXBEAMS)

!  Fourier-component solutions at ALL levels
!     ! @@ Rob Spurr, 17 July 2013, Version 2.2 --> Optional Output at ALL LEVELS

      REAL(kind=dp), INTENT(INOUT) :: RADLEVEL_F_DN (MAX_USER_STREAMS,MAXBEAMS,0:MAXLAYERS)

!  Cumulative source terms

      REAL(kind=dp), INTENT(OUT) :: CUMSOURCE_DN(MAX_USER_STREAMS,0:MAXLAYERS)

!  local variables
!  ---------------

!  Help variables

      LOGICAL       :: LUOGSS
      INTEGER       :: UM, NC, N, M, LUM
      REAL(kind=dp) :: LAYERSOURCE(MAX_USER_STREAMS)
      REAL(kind=dp) :: SHOM, SFOR2, TM

!  Local user index !@@ Observation Geometry choice, 12/21/12
!  Local user flag  !@@ Observation Geometry choice, 12/21/12

      LUM    = 1
      LUOGSS = DO_USER_OBSGEOMS .and. DO_SOLAR_SOURCES

!  Dummy (for debug)

      M = FOURIER_COMPONENT

!  Zero all Fourier components - New rule, better for safety
!    Only did this for components close to zenith (formerly)
!            !@@ 2p1, Observation Geometry choice, 12/21/12

      IF ( LUOGSS ) THEN
        INTENSITY_F_DN(LUM,IPARTIC) = 0.0d0
      ELSE
        DO UM = 1, N_USER_STREAMS
          INTENSITY_F_DN(UM,IPARTIC) = 0.0d0
        ENDDO
      ENDIF

!  Initialize recursion for user-defined stream angles only
!            !@@ 2p1, Observation Geometry choice, 12/21/12

      NC = 0
      IF ( LUOGSS ) THEN
        CUMSOURCE_DN(IPARTIC,NC) = 0.0d0
      ELSE
        DO UM = 1, N_USER_STREAMS
          CUMSOURCE_DN(UM,NC) = 0.0d0
        ENDDO
      ENDIF

!  Cumulative source terms to layer NUT (user-defined stream angles only)
!    1. Get layer source terms
!    2. Find cumulative source term

      DO N = 1, NLAYERS
        NC = N

!  Homogeneous solutions !@@ 2p1, Observation Geometry choice, 12/21/12

        IF ( LUOGSS ) THEN
          SHOM = LCON(N) * U_XNEG(IPARTIC,N) * HMULT_1(IPARTIC,N) + &
                 MCON(N) * U_XPOS(IPARTIC,N) * HMULT_2(IPARTIC,N)
          LAYERSOURCE(IPARTIC) = SHOM
        ELSE
          DO UM = 1, N_USER_STREAMS
            SHOM = LCON(N) * U_XNEG(UM,N) * HMULT_1(UM,N) + &
                   MCON(N) * U_XPOS(UM,N) * HMULT_2(UM,N)
            LAYERSOURCE(UM) = SHOM
          ENDDO
        ENDIF

!  Add thermal emission term (direct and diffuse)
!     Modulus 4.pi if solar sources are included (taken care of earlier)
!       !@@ Observation Geometry choice, 12/21/12

        IF ( DO_INCLUDE_THERMEMISS ) THEN
          TM = 1.0_dp ; IF ( DO_SOLAR_SOURCES ) TM = 1.0_dp/PI4
          IF ( LUOGSS ) THEN
            LAYERSOURCE(IPARTIC) = LAYERSOURCE(IPARTIC) + LAYER_TSUP_DN(IPARTIC,N)*TM
          ELSE
            DO UM = 1, N_USER_STREAMS
              LAYERSOURCE(UM) = LAYERSOURCE(UM) + LAYER_TSUP_DN(UM,N)*TM
            ENDDO
          ENDIF
        ENDIF

!  Add solar source term
!       !@@ Observation Geometry choice, 12/21/12

        IF ( DO_SOLAR_SOURCES ) THEN
          IF ( DO_USER_OBSGEOMS ) THEN
            SFOR2 =  EMULT_DN(LUM,N,IPARTIC) * U_WNEG2(IPARTIC,N)
            LAYERSOURCE(IPARTIC) = LAYERSOURCE(IPARTIC) + SFOR2
          ELSE
            DO UM = 1, N_USER_STREAMS
              SFOR2 =  EMULT_DN(UM,N,IPARTIC) * U_WNEG2(UM,N)
              LAYERSOURCE(UM) = LAYERSOURCE(UM) + SFOR2
            ENDDO
          ENDIF
        ENDIF

!   This is not required-----------------------No SS correction
!        IF ( .NOT.DO_MSMODE_LIDORT ) THEN
!          DO UM = 1, N_USER_STREAMS
!            SFOR1 = U_WNEG1(UM,N) * EMULT_DN(UM,N,IPARTIC)
!            LAYERSOURCE(UM) = LAYERSOURCE(UM) + SFOR1
!          ENDDO
!        ENDIF

!  Upgrade source function. 
!        !@@ 2p1, Observation Geometry choice, 12/21/12
!        !@@ 2p2, Set All-level outputs, 7/17/13

        IF ( LUOGSS ) THEN
          CUMSOURCE_DN(IPARTIC,NC) = LAYERSOURCE(IPARTIC) + &
                T_DELT_USERM(N,IPARTIC)*CUMSOURCE_DN(IPARTIC,NC-1)
          if (DO_2S_LEVELOUT)RADLEVEL_F_DN(LUM,IPARTIC,N) = FLUX_MULTIPLIER * CUMSOURCE_DN(IPARTIC,NC)
        ELSE
          DO UM = 1, N_USER_STREAMS
            CUMSOURCE_DN(UM,NC) = LAYERSOURCE(UM) + &
                         T_DELT_USERM(N,UM)*CUMSOURCE_DN(UM,NC-1)
            if (DO_2S_LEVELOUT)RADLEVEL_F_DN(UM,IPARTIC,N) = FLUX_MULTIPLIER * CUMSOURCE_DN(UM,NC)
          ENDDO
        ENDIF

      ENDDO

!    !@@ 2p1, Observation Geometry choice, 12/21/12

      IF ( LUOGSS ) THEN
        INTENSITY_F_DN(LUM,IPARTIC) = FLUX_MULTIPLIER * CUMSOURCE_DN(IPARTIC,NC)
      ELSE
        DO UM = 1, N_USER_STREAMS
          INTENSITY_F_DN(UM,IPARTIC) = FLUX_MULTIPLIER * CUMSOURCE_DN(UM,NC)
        ENDDO
      ENDIF

!  Finish

      RETURN
END SUBROUTINE TWOSTREAM_DNUSER_INTENSITY

!

SUBROUTINE TWOSTREAM_FLUXES &
            ( MAXBEAMS, MAXLAYERS, MAXTHREADS, DO_UPWELLING, DO_DNWELLING, & ! Input  Dimensions, flags
              IBEAM, NLAYERS, THREAD, PI4, STREAM_VALUE, FLUX_MULTIPLIER,  & ! Input Control
              LCON_XVEC, MCON_XVEC, EIGENTRANS, WUPPER, WLOWER,            & ! Input 2-stream solution
              FLUXES_TOA, FLUXES_BOA )                                       ! Output

!  New routine 11/5/13. Diffuse Fluxes at TOA and BOA

      implicit none

!  Precision

      INTEGER, PARAMETER :: dp = KIND( 1.0D0 )

!  Input variables
!  ---------------

!  Dimensions (2p2, add MAXLAYERS)

      INTEGER, INTENT(IN)        :: MAXBEAMS, MAXTHREADS, MAXLAYERS

!  Control

      LOGICAL, INTENT(IN)        :: DO_UPWELLING, DO_DNWELLING

!  Thread, beam, nlayers

      INTEGER, INTENT(IN)        :: NLAYERS, THREAD, IBEAM

!  multiplier, 4pi

      REAL(kind=dp), INTENT(IN)  :: FLUX_MULTIPLIER, PI4

!  Stream value

      REAL(kind=dp), INTENT(IN)  :: STREAM_VALUE

!  Eigen-Transmittance

      REAL(kind=dp), INTENT(IN) :: EIGENTRANS(MAXLAYERS)

!  Solution constants of integration multiplied by eigensolutions

      REAL(kind=dp), INTENT(IN) :: LCON_XVEC(2,MAXLAYERS)
      REAL(kind=dp), INTENT(IN) :: MCON_XVEC(2,MAXLAYERS)

!  Solutions at layer boundaries

      REAL(kind=dp), INTENT(IN) :: WUPPER(2,MAXLAYERS)
      REAL(kind=dp), INTENT(IN) :: WLOWER(2,MAXLAYERS)

!  Flux output (already initialized here)
!     ! @@ Rob Spurr, 05 November 2013, Version 2.3 --> Flux Output

      REAL(kind=dp), INTENT(INOUT) :: FLUXES_TOA(MAXBEAMS,2,MAXTHREADS)
      REAL(kind=dp), INTENT(INOUT) :: FLUXES_BOA(MAXBEAMS,2,MAXTHREADS)

!  Local variables

      INTEGER :: N
      REAL(kind=dp) :: PI2, SHOM, SPAR, QUADINTENS_TOA, QUADINTENS_BOA

!  upwelling Flux at TOA

      PI2 = 0.5d0 * PI4
      if ( DO_UPWELLING ) THEN
         N = 1
         SHOM = LCON_XVEC(2,N) + MCON_XVEC(2,N) * EIGENTRANS(N)
         SPAR = WUPPER(2,N)
         QUADINTENS_TOA = FLUX_MULTIPLIER * ( SPAR + SHOM )
         FLUXES_TOA(IBEAM,1,THREAD) = 0.5d0 * QUADINTENS_TOA
         FLUXES_TOA(IBEAM,2,THREAD) = PI2 * STREAM_VALUE * QUADINTENS_TOA
      endif

!  Downwelling Flux at BOA

      if ( DO_DNWELLING ) THEN
         N = NLAYERS
         SHOM = LCON_XVEC(1,N) * EIGENTRANS(N) + MCON_XVEC(1,N)
         SPAR = WLOWER(1,N)
         QUADINTENS_BOA = FLUX_MULTIPLIER * ( SPAR + SHOM )
         FLUXES_BOA(IBEAM,1,THREAD) = 0.5d0 * QUADINTENS_BOA
         FLUXES_BOA(IBEAM,2,THREAD) = PI2 * STREAM_VALUE * QUADINTENS_BOA
      endif

!  Finish

      RETURN
END SUBROUTINE TWOSTREAM_FLUXES

!

SUBROUTINE TWOSTREAM_CONVERGE &
       ( MAXBEAMS, MAX_USER_STREAMS, MAX_USER_RELAZMS,   & ! Dimensions
         MAX_GEOMETRIES, MAXTHREADS, MAXLAYERS,          & ! Dimensions ! @@ 2p2
         DO_UPWELLING, DO_DNWELLING, DO_2S_LEVELOUT,     & ! Inputs     ! @@ 2p2
         NLAYERS, THREAD, IBEAM, FOURIER_COMPONENT,      & ! Inputs     ! @@ 2p2
         N_USER_STREAMS, N_USER_RELAZMS, AZMFAC, UMOFF,  & ! Inputs
         INTENSITY_F_UP,  INTENSITY_F_DN,                & ! Inputs
         RADLEVEL_F_UP,   RADLEVEL_F_DN,                 & ! Inputs     ! @@ 2p2
         INTENSITY_TOA, INTENSITY_BOA,                   & ! In/Out
         RADLEVEL_UP,   RADLEVEL_DN   )                    ! In/Out     ! @@ 2p2

!  Alterations for version 2.2, 17 July 2013

      implicit none

!  Precision

      INTEGER, PARAMETER :: dp = KIND( 1.0D0 )

!  Input variables
!  ---------------

!  Dimensions (2p2, add MAXLAYERS)

      INTEGER, INTENT(IN)        :: MAXBEAMS, MAX_USER_STREAMS, MAX_USER_RELAZMS
      INTEGER, INTENT(IN)        :: MAX_GEOMETRIES, MAXTHREADS, MAXLAYERS

!  Control
!     ! @@ Rob Spurr, 17 July 2013, Version 2.2, Levelout flag

      LOGICAL, INTENT(IN)        :: DO_UPWELLING, DO_DNWELLING
      LOGICAL, INTENT(IN)        :: DO_2S_LEVELOUT

!  SS control, not required in this streamlined version
!      LOGICAL, INTENT(IN) :: DO_SSFULL, DO_SSCORR_OUTGOING, DO_SSCORR_NADIR

!  Numbers

      INTEGER, INTENT(IN)        :: N_USER_STREAMS, N_USER_RELAZMS

!  Fourier component and thread, beam, nlayers (2p2, added)

      INTEGER, INTENT(IN)        :: FOURIER_COMPONENT, NLAYERS, THREAD, IBEAM

!  Local  azimuth factors

      INTEGER, INTENT(IN)        :: UMOFF ( MAXBEAMS, MAX_USER_STREAMS )
      REAL(kind=dp), INTENT(IN)  :: AZMFAC(MAX_USER_STREAMS,MAXBEAMS,MAX_USER_RELAZMS)

!  User-defined solutions

      REAL(kind=dp), INTENT(IN)  :: INTENSITY_F_UP(MAX_USER_STREAMS,MAXBEAMS)
      REAL(kind=dp), INTENT(IN)  :: INTENSITY_F_DN(MAX_USER_STREAMS,MAXBEAMS)

!  Fourier-component solutions at ALL levels
!     ! @@ Rob Spurr, 17 July 2013, Version 2.2 --> Optional Output at ALL LEVELS

      REAL(kind=dp), INTENT(IN) :: RADLEVEL_F_UP (MAX_USER_STREAMS,MAXBEAMS,0:MAXLAYERS)
      REAL(kind=dp), INTENT(IN) :: RADLEVEL_F_DN (MAX_USER_STREAMS,MAXBEAMS,0:MAXLAYERS)

!  Single scatter solutions, Not required here
!      REAL(kind=dp), INTENT(IN)  :: INTENSITY_SS_UP(MAX_GEOMETRIES)
!      REAL(kind=dp), INTENT(IN)  :: INTENSITY_SS_DN(MAX_GEOMETRIES)

!  Output
!  ------

!  TOA and BOA output

      REAL(kind=dp), INTENT(INOUT) :: INTENSITY_TOA(MAX_GEOMETRIES,MAXTHREADS)
      REAL(kind=dp), INTENT(INOUT) :: INTENSITY_BOA(MAX_GEOMETRIES,MAXTHREADS)

!  output solutions at ALL levels
!     ! @@ Rob Spurr, 17 July 2013, Version 2.2 --> Optional Output at ALL LEVELS

      REAL(kind=dp), INTENT(INOUT) :: RADLEVEL_UP (MAX_GEOMETRIES,0:MAXLAYERS,MAXTHREADS)
      REAL(kind=dp), INTENT(INOUT) :: RADLEVEL_DN (MAX_GEOMETRIES,0:MAXLAYERS,MAXTHREADS)

!  Local variables
!  ---------------

      INTEGER       :: I, UA, V
      REAL(kind=dp) :: TOLD, TAZM

!  ###################
!  Fourier 0 component
!  ###################

      IF ( FOURIER_COMPONENT.EQ.0 ) THEN

!  Copy DIFFUSE Fourier component at all output angles and optical depths
!    If no SSCORR and no DBCORR, then two options apply:
!     (a) Convergence on RADIANCE = DIFFUSE + SSTRUNCATED + DBTRUNCATED
!              (full radiance, no SS correction, no DB correction)
!     (b) Convergence on RADIANCE = DIFFUSE alone (MS only mode)
!              (SSTRUNCATED + DBTRUNCATED do not get calculated)

!  Code only for the NON-SSFULL case
!        IF ( .not. DO_SSFULL ) THEN

!     ! @@ Rob Spurr, 17 July 2013, Version 2.2 --> Optional Output at ALL LEVELS

          DO I = 1, N_USER_STREAMS
            DO UA = 1, N_USER_RELAZMS
              V = UMOFF(IBEAM,I) + UA
              IF ( DO_UPWELLING ) THEN
                INTENSITY_TOA(V,THREAD) = INTENSITY_F_UP(I,IBEAM)
                IF ( DO_2S_LEVELOUT ) RADLEVEL_UP(V,0:NLAYERS,THREAD) = RADLEVEL_F_UP(I,IBEAM,0:NLAYERS)
              ENDIF
              IF ( DO_DNWELLING ) THEN
                INTENSITY_BOA(V,THREAD) = INTENSITY_F_DN(I,IBEAM)
                IF ( DO_2S_LEVELOUT ) RADLEVEL_DN(V,0:NLAYERS,THREAD) = RADLEVEL_F_DN(I,IBEAM,0:NLAYERS)
              ENDIF
            ENDDO
          ENDDO

!  Commented out in the streamlined version - NO SS OUTPUT
!        IF ( DO_SSFULL ) THEN
!           DO I = 1, N_USER_STREAMS
!            DO UA = 1, N_USER_RELAZMS
!              V = UMOFF(IBEAM,I) + UA
!              IF ( DO_UPWELLING ) THEN
!                INTENSITY_TOA(V,THREAD) = 0.0d0
!              ENDIF
!              IF ( DO_DNWELLING ) THEN
!                INTENSITY_BOA(V,THREAD) = 0.0d0
!              ENDIF
!            ENDDO
!          ENDDO
!        ENDIF

!    Add the single scatter component if flagged
!  Commented out in the streamlined version - NO SS OUTPUT HERE
!        IF ( DO_SSFULL.OR.DO_SSCORR_NADIR.OR.DO_SSCORR_OUTGOING ) THEN
!           DO I = 1, N_USER_STREAMS
!            DO UA = 1, N_USER_RELAZMS
!              V = UMOFF(IBEAM,I) + UA
!              IF ( DO_UPWELLING ) THEN
!                INTENSITY_TOA(V,THREAD) = 
!     &             INTENSITY_TOA(V,THREAD) + INTENSITY_SS_UP(V)
!              ENDIF
!              IF ( DO_DNWELLING ) THEN
!                INTENSITY_BOA(V,THREAD) = 
!     &             INTENSITY_BOA(V,THREAD) + INTENSITY_SS_DN(V)
!              ENDIF
!            ENDDO
!          ENDDO
!      ENDIF

!  ######################
!  Fourier component = 1
!  ######################

      ELSE

!  No examination of convergence
!  -----------------------------

!     ! @@ Rob Spurr, 17 July 2013, Version 2.2 --> Optional Output at ALL LEVELS

        DO UA = 1, N_USER_RELAZMS
          DO I = 1, N_USER_STREAMS
            V = UMOFF(IBEAM,I) + UA
            IF ( DO_UPWELLING ) THEN
              TOLD = INTENSITY_TOA(V,THREAD)
              TAZM = AZMFAC(I,IBEAM,UA)*INTENSITY_F_UP(I,IBEAM)
              INTENSITY_TOA(V,THREAD) = TOLD + TAZM
              IF ( DO_2S_LEVELOUT ) RADLEVEL_UP(V,0:NLAYERS,THREAD) = &
                   RADLEVEL_UP(V,0:NLAYERS,THREAD) + AZMFAC(I,IBEAM,UA) * RADLEVEL_F_UP(I,IBEAM,0:NLAYERS)
            ENDIF
            IF ( DO_DNWELLING ) THEN
              TOLD = INTENSITY_BOA(V,THREAD)
              TAZM = AZMFAC(I,IBEAM,UA)*INTENSITY_F_DN(I,IBEAM)
              INTENSITY_BOA(V,THREAD) = TOLD + TAZM
              IF ( DO_2S_LEVELOUT ) RADLEVEL_DN(V,0:NLAYERS,THREAD) = &
                   RADLEVEL_DN(V,0:NLAYERS,THREAD) + AZMFAC(I,IBEAM,UA) * RADLEVEL_F_DN(I,IBEAM,0:NLAYERS)
            ENDIF
          ENDDO
        ENDDO

!  Finish Fourier

      ENDIF

!  Finish

      RETURN
END SUBROUTINE TWOSTREAM_CONVERGE

!

SUBROUTINE TWOSTREAM_CONVERGE_OBSGEO &
       ( MAXBEAMS, MAX_USER_STREAMS, MAX_USER_RELAZMS,      & ! Dimensions
         MAX_GEOMETRIES, MAXTHREADS, MAXLAYERS,             & ! Dimensions ! @@ 2p2
         DO_UPWELLING, DO_DNWELLING, DO_2S_LEVELOUT,        & ! Inputs     ! @@ 2p2
         NLAYERS, THREAD, IBEAM, FOURIER_COMPONENT, AZMFAC, & ! Inputs     ! @@ 2p2
         INTENSITY_F_UP,  INTENSITY_F_DN,                   & ! Inputs
         RADLEVEL_F_UP,   RADLEVEL_F_DN,                    & ! Inputs     ! @@ 2p2
         INTENSITY_TOA,   INTENSITY_BOA,                    & ! In/Out
         RADLEVEL_UP,     RADLEVEL_DN   )                     ! In/Out     ! @@ 2p2

!  Alterations for version 2.2, 17 July 2013

      implicit none

!  Precision

      INTEGER, PARAMETER :: dp = KIND( 1.0D0 )

!  Input variables
!  ---------------

!  Dimensions (2p2, add MAXLAYERS)

      INTEGER, INTENT(IN)        :: MAXBEAMS, MAX_USER_STREAMS, MAX_USER_RELAZMS
      INTEGER, INTENT(IN)        :: MAX_GEOMETRIES, MAXTHREADS, MAXLAYERS

!  Control
!     ! @@ Rob Spurr, 17 July 2013, Version 2.2, Levelout flag

      LOGICAL, INTENT(IN)        :: DO_UPWELLING, DO_DNWELLING
      LOGICAL, INTENT(IN)        :: DO_2S_LEVELOUT

!  SS control, not required in this streamlined version
!      LOGICAL, INTENT(IN) :: DO_SSFULL, DO_SSCORR_OUTGOING, DO_SSCORR_NADIR

!  Fourier component and thread, beam, nlayers (2p2, added)

      INTEGER, INTENT(IN)        :: FOURIER_COMPONENT, NLAYERS, THREAD, IBEAM

!  Local  azimuth factors

      REAL(kind=dp), INTENT(IN)  :: AZMFAC(MAX_USER_STREAMS,MAXBEAMS,MAX_USER_RELAZMS)

!  User-defined solutions

      REAL(kind=dp), INTENT(IN)  :: INTENSITY_F_UP(MAX_USER_STREAMS,MAXBEAMS)
      REAL(kind=dp), INTENT(IN)  :: INTENSITY_F_DN(MAX_USER_STREAMS,MAXBEAMS)

!  Fourier-component solutions at ALL levels
!     ! @@ Rob Spurr, 17 July 2013, Version 2.2 --> Optional Output at ALL LEVELS

      REAL(kind=dp), INTENT(IN) :: RADLEVEL_F_UP (MAX_USER_STREAMS,MAXBEAMS,0:MAXLAYERS)
      REAL(kind=dp), INTENT(IN) :: RADLEVEL_F_DN (MAX_USER_STREAMS,MAXBEAMS,0:MAXLAYERS)

!  Single scatter solutions, Not required here
!      REAL(kind=dp), INTENT(IN)  :: INTENSITY_SS_UP(MAX_GEOMETRIES)
!      REAL(kind=dp), INTENT(IN)  :: INTENSITY_SS_DN(MAX_GEOMETRIES)

!  Output
!  ------

!  TOA and BOA output

      REAL(kind=dp), INTENT(INOUT) :: INTENSITY_TOA(MAX_GEOMETRIES,MAXTHREADS)
      REAL(kind=dp), INTENT(INOUT) :: INTENSITY_BOA(MAX_GEOMETRIES,MAXTHREADS)

!  output solutions at ALL levels
!     ! @@ Rob Spurr, 17 July 2013, Version 2.2 --> Optional Output at ALL LEVELS

      REAL(kind=dp), INTENT(INOUT) :: RADLEVEL_UP (MAX_GEOMETRIES,0:MAXLAYERS,MAXTHREADS)
      REAL(kind=dp), INTENT(INOUT) :: RADLEVEL_DN (MAX_GEOMETRIES,0:MAXLAYERS,MAXTHREADS)

!  Local variables
!  ---------------

      INTEGER       :: LUM, LUA
      REAL(kind=dp) :: TOLD, TAZM

!  Local user indices

      LUM = 1
      LUA = 1

!  Fourier 0 component
!     ! @@ Rob Spurr, 17 July 2013, Version 2.2 --> Optional Output at ALL LEVELS

      IF ( FOURIER_COMPONENT.EQ.0 ) THEN

        IF ( DO_UPWELLING ) THEN
          INTENSITY_TOA(IBEAM,THREAD) = INTENSITY_F_UP(LUM,IBEAM)
          IF ( DO_2S_LEVELOUT ) RADLEVEL_UP(IBEAM,0:NLAYERS,THREAD) = RADLEVEL_F_UP(LUM,IBEAM,0:NLAYERS)
        ENDIF
        IF ( DO_DNWELLING ) THEN
          INTENSITY_BOA(IBEAM,THREAD) = INTENSITY_F_DN(LUM,IBEAM)
          IF ( DO_2S_LEVELOUT ) RADLEVEL_DN(IBEAM,0:NLAYERS,THREAD) = RADLEVEL_F_DN(LUM,IBEAM,0:NLAYERS)
        ENDIF

!  Fourier component = 1
!     ! @@ Rob Spurr, 17 July 2013, Version 2.2 --> Optional Output at ALL LEVELS

      ELSE
        IF ( DO_UPWELLING ) THEN
          TOLD = INTENSITY_TOA(IBEAM,THREAD)
          TAZM = AZMFAC(LUM,IBEAM,LUA)*INTENSITY_F_UP(LUM,IBEAM)
          INTENSITY_TOA(IBEAM,THREAD) = TOLD + TAZM
          IF ( DO_2S_LEVELOUT ) RADLEVEL_UP(IBEAM,0:NLAYERS,THREAD) = &
                   RADLEVEL_UP(IBEAM,0:NLAYERS,THREAD) + AZMFAC(LUM,IBEAM,LUA) * RADLEVEL_F_UP(LUM,IBEAM,0:NLAYERS)
        ENDIF
        IF ( DO_DNWELLING ) THEN
          TOLD = INTENSITY_BOA(IBEAM,THREAD)
          TAZM = AZMFAC(LUM,IBEAM,LUA)*INTENSITY_F_DN(LUM,IBEAM)
          INTENSITY_BOA(IBEAM,THREAD) = TOLD + TAZM
          IF ( DO_2S_LEVELOUT ) RADLEVEL_DN(IBEAM,0:NLAYERS,THREAD) = &
                   RADLEVEL_DN(IBEAM,0:NLAYERS,THREAD) + AZMFAC(LUM,IBEAM,LUA) * RADLEVEL_F_DN(LUM,IBEAM,0:NLAYERS)
        ENDIF

!  Finish Fourier

      ENDIF

!  Finish

      RETURN
END SUBROUTINE TWOSTREAM_CONVERGE_OBSGEO

end module twostream_intensity_m

