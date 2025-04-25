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
! #   setups                                                    #
! #                                                             #
! #          TWOSTREAM_THERMALSETUP                             #
! #                                                             #
! #   discrete ordinate particular integral                     #
! #                                                             #
! #          TWOSTREAM_THERMAL_CLSOLUTION                       #
! #          TWOSTREAM_THERMAL_GFSOLUTION                       #
! #                                                             #
! #   postprocessing source terms                               #
! #                                                             #
! #          TWOSTREAM_THERMAL_CLSTERMS                         #
! #          TWOSTREAM_THERMAL_GFSTERMS                         #
! #                                                             #
! #     Specialist solution ( 5/9/24 )                          #
! #                                                             #
! #          TWOSTREAM_THERMAL_HMSOLUTION                       #
! #                                                             #
! ###############################################################

!  Rob Fix 5/24/17 for 2p4a, reintroduce classical PI subroutines:

!           TWOSTREAM_THERMAL_CLSOLUTION
!           TWOSTREAM_THERMAL_CLSTERMS

!  4/12/24. Add conservative scattering control
!    -- For situations where some (but not all) layers have conservative scattering

!  5/8/24. Add DO_HEMISPHER_MEAN option (new subroutine)

module twostream_thermalsup_m

PUBLIC

contains

SUBROUTINE TWOSTREAM_THERMALSETUP &
          ( MAXLAYERS, NLAYERS, CONSSCAT,       & ! Inputs
            OMEGA_TOTAL, DELTAU_VERT, BB_INPUT, & ! inputs
            THERMCOEFFS, TCOM1 )                  ! Outputs

      implicit none

!  precision

      INTEGER, PARAMETER :: dp     = KIND( 1.0D0 )

!  subroutine arguments
!  --------------------

!  Dimensions

      INTEGER, INTENT(IN)        :: MAXLAYERS

!  Number of layers

      INTEGER, INTENT(IN)        :: NLAYERS

!  4/12/24. Add conservative scattering control

      LOGICAL, intent(in)        :: CONSSCAT(MAXLAYERS)

!  Optical properties

      REAL(kind=dp), INTENT (IN) :: OMEGA_TOTAL ( MAXLAYERS )
      REAL(kind=dp), INTENT (IN) :: DELTAU_VERT ( MAXLAYERS )

!  Thermal input

      REAL(kind=dp), INTENT (IN) :: BB_INPUT    ( 0:MAXLAYERS )

!  Output Help variables

      REAL (kind=dp), INTENT (OUT) :: TCOM1 ( MAXLAYERS, 2 )
      REAL (kind=dp), INTENT (OUT) :: THERMCOEFFS  ( MAXLAYERS, 2 )

!  Local variables

      INTEGER        :: N
      REAL (kind=dp) :: OMEGAS1

!  thermal coefficients and powers

      DO N = 1, NLAYERS
        THERMCOEFFS(N,1)  = BB_INPUT(N-1)
        THERMCOEFFS(N,2)  = (BB_INPUT(N)-BB_INPUT(N-1))/DELTAU_VERT(N)
      END DO

!  Auxiliary output

      DO N = 1, NLAYERS
        IF ( CONSSCAT(N) ) then
          TCOM1(N,1:2) = 0.0_dp
        ELSE
          OMEGAS1    = 1.0_dp - OMEGA_TOTAL(N)
          TCOM1(N,1) = THERMCOEFFS(N,1) * OMEGAS1
          TCOM1(N,2) = THERMCOEFFS(N,2) * OMEGAS1
        ENDIF
      END DO

!  Finish

      RETURN
END SUBROUTINE TWOSTREAM_THERMALSETUP

!

SUBROUTINE TWOSTREAM_THERMAL_HMSOLUTION &
      ( MAXLAYERS, NLAYERS, CONSSCAT, DELTAU_VERT, & ! Inputs
        THERMCOEFFS, GAMMA1, GAMMA2,               & ! Inputs
        T_WUPPER, T_WLOWER )                         ! Outputs

!  5/8/24. New subroutine for the DO_HEMISPHER_MEAN case

      implicit none

!  precision and parameters

      INTEGER      , PARAMETER :: dp   = KIND( 1.0D0 )
      REAL(kind=dp), parameter :: ZERO = 0.0_dp, ONE = 1.0_dp

!  subroutine arguments
!  --------------------

!  Numbers

      INTEGER, INTENT(IN)        :: MAXLAYERS, NLAYERS

!  deltau

      REAL (kind=dp), INTENT (IN) :: DELTAU_VERT ( MAXLAYERS )

!  Thermal coefficients

      REAL (kind=dp), INTENT (IN) :: THERMCOEFFS ( MAXLAYERS, 2 )

!  4/12/24. Add conservative scattering control

      LOGICAL, intent(in)         :: CONSSCAT(MAXLAYERS)

!  local matrices from eigenvalue computation
!      REAL(kind=dp), INTENT(IN)   :: DAB(MAXLAYERS)

!  local matrices from Hemisphere-mean  computation

      REAL(kind=dp), INTENT(IN)   :: GAMMA1(MAXLAYERS)
      REAL(kind=dp), INTENT(IN)   :: GAMMA2(MAXLAYERS)

!  Output variables
!  ----------------

!  Thermal solution at layer boundaries

      REAL(kind=dp) , INTENT (OUT) ::  T_WUPPER ( 2, MAXLAYERS )
      REAL(kind=dp) , INTENT (OUT) ::  T_WLOWER ( 2, MAXLAYERS )

!  Local variables
! ---------------

      INTEGER       :: N
      REAL(kind=dp) :: DELFAC, TMAT, HVEC1, JVEC1

!  ZERO THE BOUNDARY LAYER VALUES

      T_WUPPER(1:2,1:NLAYERS) = zero
      T_WLOWER(1:2,1:NLAYERS) = zero

! ---------------------------------------
! CLASSICAL SOLUTIONS FOR ALL LAYERS
! ---------------------------------------

      DO N = 1, NLAYERS

!  4/12/24. Add conservative scattering control

        IF ( CONSSCAT(N) ) CYCLE

!  2-stream Solutions (hemispheric mean)

!        TMAT   = - DAB(N)
        TMAT   = GAMMA1(n) + GAMMA2(n)
        HVEC1  = THERMCOEFFS(N,1)
        JVEC1  = THERMCOEFFS(N,2) / TMAT
        DELFAC = THERMCOEFFS(N,2) * DELTAU_VERT(N)
 
!write(*,*)n,gamma1(1),gamma2(n),tmat

!  Stream Solution, 1 = down, 2 = up. Values at layer boundaries

        T_WUPPER(1,N) = HVEC1 - JVEC1
        T_WUPPER(2,N) = HVEC1 + JVEC1
        T_WLOWER(1,N) = T_WUPPER(1,N) + DELFAC
        T_WLOWER(2,N) = T_WUPPER(2,N) + DELFAC

!  END LAYER LOOP

      ENDDO

!  finish

      RETURN
END SUBROUTINE TWOSTREAM_THERMAL_HMSOLUTION

!

SUBROUTINE TWOSTREAM_THERMAL_CLSOLUTION  &
      ( MAXLAYERS, MAX_USER_STREAMS, DO_UPWELLING, DO_DNWELLING,       & ! Inputs
        CONSSCAT, NLAYERS, N_USER_STREAMS, STREAM_VALUE, USER_STREAMS, & ! Inputs
        OMEGA, ASYMM, SAB, DAB, DELTAU_VERT, THERMCOEFFS,              & ! Inputs
        T_WUPPER, T_WLOWER, U_TPOS2, U_TNEG2 )                           ! Outputs

      implicit none

!  precision and parameters

      INTEGER      , PARAMETER :: dp   = KIND( 1.0D0 )
      REAL(kind=dp), parameter :: ZERO = 0.0_dp, ONE = 1.0_dp

!  subroutine arguments
!  --------------------

!  Dimensions

      INTEGER, INTENT(IN)        :: MAXLAYERS
      INTEGER, INTENT(IN)        :: MAX_USER_STREAMS

!  Flags

      LOGICAL, INTENT (IN)        :: DO_UPWELLING
      LOGICAL, INTENT (IN)        :: DO_DNWELLING

!  4/12/24. Add conservative scattering control

      LOGICAL, intent(in)         :: CONSSCAT(MAXLAYERS)

!  Numbers

      INTEGER, INTENT (IN)        :: NLAYERS
      INTEGER, INTENT (IN)        :: N_USER_STREAMS

!  Stream values

      REAL(kind=dp), INTENT(IN)   :: STREAM_VALUE
      REAL(kind=dp), INTENT(IN)   :: USER_STREAMS  ( MAX_USER_STREAMS )

!  OMEGA and ASYMM

      REAL(kind=dp), INTENT(IN)   :: OMEGA ( MAXLAYERS )
      REAL(kind=dp), INTENT(IN)   :: ASYMM ( MAXLAYERS )

!  local matrices from eigenvalue computation

      REAL(kind=dp), INTENT(IN)   :: SAB(MAXLAYERS), DAB(MAXLAYERS)

!  deltau

      REAL (kind=dp), INTENT (IN) :: DELTAU_VERT ( MAXLAYERS )

!  Thermal coefficients

      REAL (kind=dp), INTENT (IN) :: THERMCOEFFS ( MAXLAYERS, 2 )

!  Output variables
!  ----------------

!  Thermal solution at layer boundaries

      REAL(kind=dp) , INTENT (OUT) ::  T_WUPPER ( 2, MAXLAYERS )
      REAL(kind=dp) , INTENT (OUT) ::  T_WLOWER ( 2, MAXLAYERS )

!  User solution

      REAL(kind=dp) , INTENT (OUT) ::  U_TPOS2 ( MAX_USER_STREAMS, MAXLAYERS, 2 )
      REAL(kind=dp) , INTENT (OUT) ::  U_TNEG2 ( MAX_USER_STREAMS, MAXLAYERS, 2 )

!  Local variables
! ---------------

      INTEGER       :: UM, S, N, I
      REAL(kind=dp) :: TVEC1(2), TVEC2(2), TCOM1(2), HELP
      REAL(kind=dp) :: POS2, NEG2, TERM1, TERM2, TMAT, HVEC1, HVEC2, JVEC1
      REAL(kind=dp) :: WT_HELP(0:1,2), T_HELP2(0:1,2), HMU_STREAM, OMEGA_MOM

!  ZERO THE BOUNDARY LAYER VALUES

      T_WUPPER(1:2,1:NLAYERS) = zero
      T_WLOWER(1:2,1:NLAYERS) = zero

!  4/12/24. ZERO the  other outputs

      U_TPOS2(:,1:NLAYERS,:) = zero
      U_TNEG2(:,1:NLAYERS,:) = zero

! ---------------------------------------
! CLASSICAL SOLUTIONS FOR ALL LAYERS
! ---------------------------------------

      DO N = 1, NLAYERS

!  4/12/24. Add conservative scattering control

        IF ( CONSSCAT(N) ) CYCLE

!  2-stream Solutions

        HELP = one - OMEGA(N)
        TCOM1(1:2) = THERMCOEFFS(N,1:2) * help
        TERM1 = 2.0_dp * TCOM1(1)
        TERM2 = 2.0_dp * TCOM1(2)
        TMAT  = SAB(N) * STREAM_VALUE
        HVEC1 = - TERM1 / TMAT
        HVEC2 = - TERM2 / TMAT
        TMAT  = - DAB(N)
        JVEC1 = - HVEC2 / TMAT

!  Stream Solution, 1 = down, 2 = up. Values at layer boundaries

        TVEC1(1) = 0.5_dp * (HVEC1 + JVEC1)
        TVEC1(2) = 0.5_dp * (HVEC1 - JVEC1)
        TVEC2(1) = 0.5_dp * HVEC2
        TVEC2(2) = TVEC2(1)
        DO I = 1, 2
          T_WUPPER(I,N) = TVEC1(I)
          T_WLOWER(I,N) = TVEC1(I) + TVEC2(I) * DELTAU_VERT(N)
        ENDDO

!  Checks out with Toon Code in Full-quadrature mode
!write(*,*)'Top/bottom of layer', T_WUPPER(1:2,n), T_WLOWER(1:2,n)
 
!  USER SOLUTIONS

        OMEGA_MOM = 3.0d0 * OMEGA(N) * ASYMM(N)
        HMU_STREAM = STREAM_VALUE * 0.5d0

        WT_help(0,1) = ( TVEC1(2) + TVEC1(1) ) * 0.5d0
        WT_help(1,1) = ( TVEC1(2) - TVEC1(1) ) * HMU_STREAM
        WT_help(0,2) = ( TVEC2(2) + TVEC2(1) ) * 0.5d0
        WT_help(1,2) = ( TVEC2(2) - TVEC2(1) ) * HMU_STREAM
        T_help2(0,1:2)  = wt_help(0,1:2) * omega(n)
        T_help2(1,1:2)  = wt_help(1,1:2) * omega_mom

        IF ( DO_UPWELLING ) THEN
          DO UM = 1, N_USER_STREAMS
            DO S = 1, 2
              pos2 = T_help2(0,s) + T_help2(1,s)* user_streams(um)
              U_TPOS2(UM,N,S) = pos2
            ENDDO
          ENDDO
        ENDIF

        IF ( DO_DNWELLING ) THEN
          DO UM = 1, N_USER_STREAMS
            DO S = 1, 2
              neg2 = T_help2(0,s) - T_help2(1,s)* user_streams(um)
              U_TNEG2(UM,N,S) = neg2
            ENDDO
          ENDDO
        ENDIF

!write(*,*)n,t_wlower(1:2,n)

!  END LAYER LOOP

      ENDDO

!  finish

      RETURN
END SUBROUTINE TWOSTREAM_THERMAL_CLSOLUTION

!

SUBROUTINE TWOSTREAM_THERMAL_GFSOLUTION &
      ( MAXLAYERS, NLAYERS, OMEGA, DELTAUS, THERMCOEFFS,             & ! Inputs
        CONSSCAT, TCUTOFF, EIGENVALUE, EIGENTRANS, XPOS, NORM_SAVED, & ! Input
        T_C_PLUS, T_C_MINUS, TTERM_SAVE, T_WUPPER, T_WLOWER )          ! Outputs

      implicit none

!  precision and parameters

      INTEGER      , PARAMETER :: dp   = KIND( 1.0D0 )
      REAL(kind=dp), parameter :: ZERO = 0.0_dp, ONE = 1.0_dp

!  subroutine arguments
!  --------------------

!  Dimensions

      INTEGER, INTENT(IN)        :: MAXLAYERS

!  Numbers

      INTEGER, INTENT (IN)        :: NLAYERS

!  4/12/24. Add conservative scattering control

      LOGICAL, intent(in)         :: CONSSCAT(MAXLAYERS)

!  OMEGA and DELTAUS

      REAL(kind=dp), INTENT(IN)   :: OMEGA ( MAXLAYERS )
      REAL(kind=dp), INTENT (IN)  :: DELTAUS ( MAXLAYERS )

!  Help variable

      REAL(kind=dp), INTENT (IN) :: THERMCOEFFS ( MAXLAYERS, 2 )

!  Thermal Cutoff (actually a layer optical thickness minimum)
!     Rob, introduced 14 May 2015, following 2p3 implementation (2014)
!    Solutions are avoided for optically thin layers

      REAL(kind=dp), INTENT (IN) :: TCUTOFF

!  Transmittance factors and eigenvalues

      REAL(kind=dp), intent(in)  :: EIGENVALUE (MAXLAYERS)
      REAL(kind=dp), intent(in)  :: EIGENTRANS (MAXLAYERS)

!  Eigenvector solutions

      REAL(kind=dp), intent(in)  :: XPOS (2,MAXLAYERS)

!  Saved quantities for the Green function solution

      REAL(kind=dp), intent(in)  :: NORM_SAVED (MAXLAYERS)

!  Output variables
!  ----------------

!  Saved quantities for the Green function solution

      REAL(kind=dp), intent(out) :: TTERM_SAVE (MAXLAYERS)
      REAL(kind=dp), intent(out) :: T_C_MINUS (MAXLAYERS,0:2)
      REAL(kind=dp), intent(out) :: T_C_PLUS  (MAXLAYERS,0:2)

!  Thermal solution at layer boundaries

      REAL(kind=dp) , INTENT (OUT) ::  T_WUPPER ( 2, MAXLAYERS )
      REAL(kind=dp) , INTENT (OUT) ::  T_WLOWER ( 2, MAXLAYERS )

!  Local variables
! ---------------

      INTEGER       :: N
      REAL(kind=dp) :: TK, K1, TT, t_gmult_dn, t_gmult_up
      REAL(kind=dp) :: omega1, help, sum_m, sum_p

!  ZERO THE BOUNDARY LAYER VALUES

      T_WUPPER(1:2,1:NLAYERS) = zero
      T_WLOWER(1:2,1:NLAYERS) = zero

!  4/12/24. ZERO the  other outputs

      TTERM_SAVE(1:NLAYERS)     = zero
      T_C_MINUS (1:NLAYERS,1:2) = zero
      T_C_PLUS  (1:NLAYERS,1:2) = zero

! ---------------------------------------
! Green function solutions for all layers
! ---------------------------------------

   DO n = 1, nlayers

!  4/12/24. Add conservative scattering control

      IF ( CONSSCAT(N) ) CYCLE

!  Only do this if cutoff is not in play

      if ( DELTAUS(n) .gt. TCUTOFF )  then

!  get the constant term

        omega1 = one - omega(n)
        help   = xpos(1,n)+xpos(2,n)
        tterm_save(n) = omega1 * help / norm_saved(n)
        tt = tterm_save(n)

! Green function multipliers

        k1 = one / eigenvalue(n)
        tk = eigentrans(n)
        t_c_minus(n,2)  = k1 * thermcoeffs(n,2)
        t_c_plus(n,2)   = k1 * thermcoeffs(n,2)
        t_c_minus(n,1) = k1 * ( thermcoeffs(n,1) - t_c_minus(n,2) )
        t_c_plus(n,1)  = k1 * ( thermcoeffs(n,1) + t_c_plus (n,2) )
        sum_m = t_c_minus(n,1) + t_c_minus(n,2) * deltaus(n)
        sum_p = t_c_plus (n,1) + t_c_plus(n,2)  * deltaus(n)
        t_c_minus(n,0) = - t_c_minus(n,1)
        t_c_plus(n,0)  = - sum_p
        t_gmult_dn = tt * ( tk * t_c_minus(n,0) + sum_m )
        t_gmult_up = tt * ( tk * t_c_plus (n,0) + t_c_plus(n,1) )

! Set particular integral from Green function expansion

        t_wupper(1,n) = t_gmult_up*xpos(2,n)
        t_wupper(2,n) = t_gmult_up*xpos(1,n)
        t_wlower(2,n) = t_gmult_dn*xpos(2,n)
        t_wlower(1,n) = t_gmult_dn*xpos(1,n)

!write(*,'(i2,1p6e18.8)')n,t_wupper(1:2,n),t_wlower(1:2,n),norm_saved(n),eigentrans(n)

!  End active layer condition

      ENDIF

!write(*,*)n,t_wlower(1:2,n),omega1, norm_saved(N),eigenvalue(n)

!  End layer loop

   END DO

!  finish

   RETURN
END SUBROUTINE TWOSTREAM_THERMAL_GFSOLUTION

!

SUBROUTINE TWOSTREAM_THERMAL_CLSTERMS  &
      ( MAXLAYERS, MAX_USER_STREAMS, CONSSCAT,        & ! Dimensions/CSFlags
        DO_UPWELLING, DO_DNWELLING, DO_SOLAR_SOURCES, & ! Inputs
        NLAYERS, N_USER_STREAMS, PI4, USER_STREAMS,   & ! Inputs
        T_DELT_USERM, DELTAU_VERT, U_TPOS2, U_TNEG2,  & ! Inputs
        LAYER_TSUP_UP, LAYER_TSUP_DN  )                 ! Outputs

      implicit none

!  precision and parameters

      INTEGER      , PARAMETER :: dp   = KIND( 1.0D0 )
      REAL(kind=dp), parameter :: ZERO = 0.0_dp, ONE = 1.0_dp

!  subroutine arguments
!  --------------------

!  Dimensions

      INTEGER, INTENT(IN)         :: MAXLAYERS, MAX_USER_STREAMS

!  4/12/24. Add conservative scattering control

      LOGICAL, intent(in)         :: CONSSCAT(MAXLAYERS)

!  Flags

      LOGICAL, INTENT (IN)        :: DO_UPWELLING
      LOGICAL, INTENT (IN)        :: DO_DNWELLING
      LOGICAL, INTENT (IN)        :: DO_SOLAR_SOURCES

!  Constant

      REAL(kind=dp), INTENT (IN)  :: PI4

!  Numbers

      INTEGER, INTENT (IN)        :: NLAYERS
      INTEGER, INTENT (IN)        :: N_USER_STREAMS

!  User streams and layer transmittances

      REAL(kind=dp), INTENT (IN)  :: USER_STREAMS  ( MAX_USER_STREAMS )
      REAL(kind=dp), INTENT (IN)  :: T_DELT_USERM ( MAXLAYERS, MAX_USER_STREAMS )

!  optical depth

      REAL (kind=dp), INTENT (IN) :: DELTAU_VERT ( MAXLAYERS )

!  User solutions

      REAL(kind=dp) , INTENT (IN) :: U_TPOS2 ( MAX_USER_STREAMS, MAXLAYERS, 2 )
      REAL(kind=dp) , INTENT (IN) :: U_TNEG2 ( MAX_USER_STREAMS, MAXLAYERS, 2 )

!  Output variables
!  ----------------

      REAL(kind=dp), INTENT (OUT) :: LAYER_TSUP_UP ( MAX_USER_STREAMS, MAXLAYERS )
      REAL(kind=dp), INTENT (OUT) :: LAYER_TSUP_DN ( MAX_USER_STREAMS, MAXLAYERS )

!  Local variables
! ---------------

      INTEGER       :: UM, N
      REAL(kind=dp) :: SPAR, COSMUM, FAC, SUM
      REAL(kind=dp) :: T_MULT_UP (0:2), T_MULT_DN (0:2)

!  PARTICULAR SOLUTION LAYER SOURCE TERMS ( CLASSICAL SOLUTION )
!  --------------------------------------------------------------------

!  Initialize

      LAYER_TSUP_UP = zero
      LAYER_TSUP_DN = zero

!  INITIAL MODULUS = 4.PI IF SOLAR SOURCES ARE INCLUDED

      FAC = 1.0_dp
      IF ( DO_SOLAR_SOURCES ) FAC = PI4

!  UPWELLING and DOWNWELLING WHOLE LAYER SOURCE TERMS
!   NOTE: T_DELT_USERM(N,UM) WAS INDEXED OPPOSITELY

!  4/12/24. Add conservative scattering control
!            IF ( CONSSCAT(N) ) CYCLE

      DO UM = 1, N_USER_STREAMS
        COSMUM = USER_STREAMS(UM)
        IF ( DO_UPWELLING ) THEN
          DO N = 1, NLAYERS
            IF ( CONSSCAT(N) ) CYCLE
            T_MULT_UP(2) = U_TPOS2(UM,N,2)
            T_MULT_UP(1) = U_TPOS2(UM,N,1) + COSMUM * T_MULT_UP(2)
            SUM = T_MULT_UP(1) + T_MULT_UP(2) * DELTAU_VERT(N)
            T_MULT_UP(0)   = - SUM
            SPAR = T_MULT_UP(0) * T_DELT_USERM(N,UM) + T_MULT_UP(1)
            LAYER_TSUP_UP(UM,N) = FAC * SPAR
          ENDDO
        ENDIF

        IF ( DO_DNWELLING ) THEN
          DO N = 1, NLAYERS
            IF ( CONSSCAT(N) ) CYCLE
            T_MULT_DN(2) = U_TNEG2(UM,N,2)
            T_MULT_DN(1) = U_TNEG2(UM,N,1) - COSMUM * T_MULT_DN(2)
            T_MULT_DN(0) = - T_MULT_DN(1)
            SPAR = T_MULT_DN(0) * T_DELT_USERM(N,UM)
            SPAR = SPAR + T_MULT_DN(1)
            SPAR = SPAR + T_MULT_DN(2) * DELTAU_VERT(N)
            LAYER_TSUP_DN(UM,N) = FAC * SPAR
          ENDDO
        ENDIF

      ENDDO

!  Finish

      RETURN
END SUBROUTINE TWOSTREAM_THERMAL_CLSTERMS

!

SUBROUTINE TWOSTREAM_THERMAL_GFSTERMS &
      ( MAXLAYERS, MAX_USER_STREAMS, CONSSCAT,        & ! Dimensions/CSFlags
        DO_UPWELLING, DO_DNWELLING, DO_SOLAR_SOURCES, & ! Inputs
        NLAYERS, N_USER_STREAMS, PI4, USER_STREAMS,   & ! Inputs
        TCUTOFF, T_DELT_USERM, DELTAUS,               & ! Inputs
        U_XPOS, U_XNEG, HMULT_1, HMULT_2,             & ! Inputs
        T_C_PLUS, T_C_MINUS, TTERM_SAVE,              & ! Inputs
        LAYER_TSUP_UP, LAYER_TSUP_DN  )                 ! Outputs

      implicit none

!  precision and parameters

      INTEGER      , PARAMETER :: dp   = KIND( 1.0D0 )
      REAL(kind=dp), parameter :: ZERO = 0.0_dp, ONE = 1.0_dp

!  subroutine arguments
!  --------------------

!  Dimensions

      INTEGER, INTENT(IN)         :: MAXLAYERS, MAX_USER_STREAMS

!  4/12/24. Add conservative scattering control

      LOGICAL, intent(in)         :: CONSSCAT(MAXLAYERS)

!  Flags

      LOGICAL, INTENT (IN)        :: DO_UPWELLING
      LOGICAL, INTENT (IN)        :: DO_DNWELLING
      LOGICAL, INTENT (IN)        :: DO_SOLAR_SOURCES

!  Constant

      REAL(kind=dp), INTENT (IN)  :: PI4

!  Numbers

      INTEGER, INTENT (IN)        :: NLAYERS
      INTEGER, INTENT (IN)        :: N_USER_STREAMS

!  User streams and layer transmittances

      REAL(kind=dp), INTENT (IN)  :: USER_STREAMS  ( MAX_USER_STREAMS )
      REAL(kind=dp), INTENT (IN)  :: T_DELT_USERM ( MAXLAYERS, MAX_USER_STREAMS )

!  Thermal Cutoff (actually a layer optical thickness minimum)
!     Rob, introduced 14 May 2015, following 2p3 implementation (2014)
!    Solutions are avoided for optically thin layers

      REAL(kind=dp), INTENT (IN) :: TCUTOFF

!  Optical thickness

      REAL(kind=dp), INTENT (IN)  :: DELTAUS ( MAXLAYERS )

!  Eigenvectors defined at user-defined stream angles

      REAL(kind=dp), intent(in)   :: U_XPOS (MAX_USER_STREAMS,MAXLAYERS)
      REAL(kind=dp), intent(in)   :: U_XNEG (MAX_USER_STREAMS,MAXLAYERS)

!  Integrated homogeneous solution multipliers (global, whole layer)

      REAL(kind=dp), intent(in)  :: HMULT_1 (MAX_USER_STREAMS,MAXLAYERS)
      REAL(kind=dp), intent(in)  :: HMULT_2 (MAX_USER_STREAMS,MAXLAYERS)

!  Saved quantities for the Green function solution

      REAL(kind=dp), intent(in)  :: TTERM_SAVE (MAXLAYERS)
      REAL(kind=dp), intent(in)  :: T_C_MINUS (MAXLAYERS,0:2)
      REAL(kind=dp), intent(in)  :: T_C_PLUS  (MAXLAYERS,0:2)

!  Output variables
!  ----------------

      REAL(kind=dp), INTENT (OUT) :: LAYER_TSUP_UP ( MAX_USER_STREAMS, MAXLAYERS )
      REAL(kind=dp), INTENT (OUT) :: LAYER_TSUP_DN ( MAX_USER_STREAMS, MAXLAYERS )

!  Local variables
! ---------------

      INTEGER       :: UM, N
      REAL(kind=dp) :: SUM_M, SUM_P, SU, SD, SPAR, COSMUM, FAC
      REAL(kind=dp) :: TSGM_DU(MAXLAYERS,0:2)
      REAL(kind=dp) :: TSGM_DD(MAXLAYERS,0:2)
      REAL(kind=dp) :: TSGM_UU(MAXLAYERS,0:2)
      REAL(kind=dp) :: TSGM_UD(MAXLAYERS,0:2)

!  PARTICULAR SOLUTION LAYER SOURCE TERMS ( GREEN'S FUNCTION SOLUTION )
!  --------------------------------------------------------------------

!  Initialize

    LAYER_TSUP_UP = zero
    LAYER_TSUP_DN = zero

!  INITIAL MODULUS = 4.PI IF SOLAR SOURCES ARE INCLUDED

    FAC = one
    IF ( DO_SOLAR_SOURCES ) FAC = PI4

!  UPWELLING and DOWNWELLING WHOLE LAYER SOURCE TERMS
!   NOTE: T_DELT_USERM(N,UM) WAS INDEXED OPPOSITELY
!  Rob  Fix, 5.14.15. Only get solutions for  valid non-Cutoff layers.

!  4/12/24. Add conservative scattering control
!            IF ( CONSSCAT(N) ) CYCLE

    DO UM = 1, N_USER_STREAMS
      COSMUM = USER_STREAMS(UM)

!  Upwelling

      IF ( DO_UPWELLING ) THEN
        DO N = 1, NLAYERS
          IF ( CONSSCAT(N) ) CYCLE
          if ( deltaus(n) .gt. TCUTOFF ) then
            tsgm_uu(n,2) = t_c_plus (n,2)
            tsgm_ud(n,2) = t_c_minus(n,2)
            tsgm_uu(n,1) = t_c_plus (n,1) + cosmum * tsgm_uu(n,2)
            tsgm_ud(n,1) = t_c_minus(n,1) + cosmum * tsgm_ud(n,2)
            tsgm_uu(n,0) = - tsgm_uu(n,1) - tsgm_uu(n,2) * deltaus(n)
            tsgm_ud(n,0) = - tsgm_ud(n,1) - tsgm_ud(n,2) * deltaus(n)
            su = t_c_plus(n,0)  * hmult_1(um,n) + &
                      tsgm_uu(n,0) * t_delt_userm(n,um) + tsgm_uu(n,1)
            sd = t_c_minus(n,0) * hmult_2(um,n) + &
                      tsgm_ud(n,0) * t_delt_userm(n,um) + tsgm_ud(n,1)
            spar = tterm_save(n) * ( u_xpos(um,n)*sd + u_xneg(um,n)*su )
            layer_tsup_up(um,n) = fac * spar
          ENDIF
        ENDDO
      ENDIF

      IF ( DO_DNWELLING ) THEN
        DO N = 1, NLAYERS
          IF ( CONSSCAT(N) ) CYCLE
          if ( deltaus(n) .gt. TCUTOFF ) then
            tsgm_du(n,2) = t_c_plus(n,2)
            tsgm_dd(n,2) = t_c_minus(n,2)
            tsgm_du(n,1) = t_c_plus (n,1) - cosmum * tsgm_du(n,2)
            tsgm_dd(n,1) = t_c_minus(n,1) - cosmum * tsgm_dd(n,2)
            tsgm_du(n,0) = - tsgm_du(n,1)
            tsgm_dd(n,0) = - tsgm_dd(n,1)
            sum_p = tsgm_du(n,1)  + tsgm_du(n,2) * deltaus(n)
            sum_m = tsgm_dd(n,1)  + tsgm_dd(n,2) * deltaus(n)
            su = t_c_plus(n,0)  * hmult_2(um,n) + &
                    tsgm_du(n,0) * t_delt_userm(n,um) + sum_p
            sd = t_c_minus(n,0) * hmult_1(um,n) + &
                    tsgm_dd(n,0) * t_delt_userm(n,um) + sum_m
            spar = tterm_save(n) * (  u_xneg(um,n)*sd + u_xpos(um,n)*su )
            layer_tsup_dn(um,n) = fac * spar
          ENDIF
        ENDDO
      ENDIF

    ENDDO

!  Finish

      RETURN
END SUBROUTINE TWOSTREAM_THERMAL_GFSTERMS

!  End module

end module twostream_thermalsup_m

