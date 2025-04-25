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
! #     Mark 9  : June   2014, Inverse Pentadiagonal        #
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
! #     Homogeneous solution                                    #
! #                                                             #
! #              TWOSTREAM_HOM_SOLUTION                         #
! #              TWOSTREAM_HOM_USERSOLUTION                     #
! #              TWOSTREAM_HMULT_MASTER                         #
! #                                                             #
! #     Particular integrals                                    #
! #                                                             #
! #              TWOSTREAM_GBEAM_SOLUTION                       #
! #              TWOSTREAM_CBEAM_SOLUTION                       #
! #              TWOSTREAM_CBEAM_USERSOLUTION                   #
! #              TWOSTREAM_CONSSCAT_SOLUTION                    #
! #                                                             #
! #     Specialist solutions ( 5/9/24 )                         #
! #                                                             #
! #              OTHER2S_HOM_SOLUTION                           #
! #              EDDINGTON_BEAM_SOLUTION                        #
! #                                                             #
! ###############################################################

!  Rob Fix 5/24/17 for 2p4a, reintroduce routine for classical PI
!              TWOSTREAM_CBEAM_SOLUTION
!              TWOSTREAM_CBEAM_USERSOLUTION

!  Rob Fix 5/24/17 for 2p4a, develop solution for Conservative scattering
!              TWOSTREAM_CONSSCAT_SOLUTION

!  Rob Fix 5/9/24, develop solutions for Specialist options

module twostream_solutions_m

!    Introduced for V2p4, Mark 10

   use Twostream_Taylor_m, only : Twostream_Taylor_Series_1

PUBLIC

contains


SUBROUTINE OTHER2S_HOM_SOLUTION &
          ( MAXLAYERS, DO_HEMISPHER_MEAN, DO_EDDINGTON_FLUX, & ! Inputs
            N, OMEGA_TOTAL, ASYMM_TOTAL, DELTAU_VERT,        & ! Inputs
            GAMMA1, GAMMA2, BIG_GAMMA, EIGENVALUE, EIGENTRANS, XPOS, XNEG ) ! In/Out

!  5/9/24. Homogeneous solutions for the 2 FLux calculations using
!          DO_HEMISPHER_MEAN (thermal only) or DO_EDDINGTON_FLUX (solar only)

      implicit none

!  precision and parameters

      INTEGER      , PARAMETER :: dp   = KIND( 1.0D0 )
      REAL(kind=dp), parameter :: ZERO = 0.0_dp, ONE = 1.0_dp

!  subroutine arguments
!  --------------------

!  Dimensions

      INTEGER, INTENT(IN)         :: MAXLAYERS

!  5/8/24. Introduce flag DO_HEMISPHER_MEAN for thermal regime

      Logical, intent(in)         :: DO_HEMISPHER_MEAN

!  5/9/24. Introduce flag DO_EDDINGTON_FLUX for solar flux regime

      Logical, intent(in)         :: DO_EDDINGTON_FLUX

!  Given layer index

      INTEGER, INTENT(IN)         :: N

!  OMEGA and ASYMM

      REAL(kind=dp), INTENT(IN)   :: OMEGA_TOTAL ( MAXLAYERS )
      REAL(kind=dp), INTENT(IN)   :: ASYMM_TOTAL ( MAXLAYERS )

!  optical thickness

      REAL(kind=dp), INTENT(IN)   :: DELTAU_VERT(MAXLAYERS)

!  Solutions to the homogeneous RT equations 
!  -----------------------------------------

!  local Quantities for eigenvalue computation

      REAL(kind=dp), INTENT(INOUT)  :: GAMMA1(MAXLAYERS), GAMMA2(MAXLAYERS)

!  Eigensolutions


      REAL(kind=dp), INTENT(INOUT)  :: BIG_GAMMA (MAXLAYERS)
      REAL(kind=dp), INTENT(INOUT)  :: EIGENVALUE(MAXLAYERS)
      REAL(kind=dp), INTENT(INOUT)  :: EIGENTRANS(MAXLAYERS)

!  UP and down solutions

      REAL(kind=dp), INTENT(INOUT)  :: XPOS(2,MAXLAYERS)
      REAL(kind=dp), INTENT(INOUT)  :: XNEG(2,MAXLAYERS)

!  Local

      real(kind=dp) :: HELP
      REAL(kind=dp), parameter :: MAX_TAU_QPATH = 200.0_dp

!  Set Gamma1 and Gamma2

      IF ( DO_HEMISPHER_MEAN ) THEN
        GAMMA1(n) = 2.0_dp - OMEGA_TOTAL(n) * ( 1.0_dp + ASYMM_TOTAL(N) )
        GAMMA2(n) = OMEGA_TOTAL(n) * ( 1.0_dp - ASYMM_TOTAL(N) )
      ELSE IF ( DO_EDDINGTON_FLUX ) THEN
        HELP = 3.0_dp * ASYMM_TOTAL(N)
        GAMMA1(n) = + 0.25_dp * ( 7.0_dp - OMEGA_TOTAL(n) * ( 4.0_dp + HELP ) )
        GAMMA2(n) = - 0.25_dp * ( 1.0_dp - OMEGA_TOTAL(n) * ( 4.0_dp - HELP ) )
      ENDIF

!  Set eigenvalue

      EIGENVALUE(N) = SQRT(GAMMA1(n)*GAMMA1(n) - GAMMA2(n)*GAMMA2(n))

!  Eigentrans, defined properly. [V2p3, Mark 10]

      HELP = EIGENVALUE(N)*DELTAU_VERT(N)
      IF ( HELP .GT. MAX_TAU_QPATH ) THEN
         EIGENTRANS(N) = zero
      ELSE
         EIGENTRANS(N) = EXP(-HELP)
      ENDIF

!  Solutions.. .Symmetry and Big Gamma

      BIG_GAMMA(N) = ( GAMMA1(n) - EIGENVALUE(N) ) / GAMMA2(N)
      XPOS(1,N) = 1.0_dp
      XPOS(2,N) = BIG_GAMMA(N)
      XNEG(1,N) = XPOS(2,N)
      XNEG(2,N) = XPOS(1,N)

   return
END SUBROUTINE OTHER2S_HOM_SOLUTION

!

SUBROUTINE TWOSTREAM_HOM_SOLUTION &
          ( MAXLAYERS, N, FOURIER, STREAM_VALUE, PXSQ, & ! Inputs
            OMEGA, ASYMM, DELTAU_VERT,                 & ! Inputs
            SAB, DAB, EIGENVALUE, EIGENTRANS,          & ! In/Out
            XPOS, XNEG, ALPHA, EPSON, GAMMA, NORM_SAVED )   ! In/Out

      implicit none

!  precision and parameters

      INTEGER      , PARAMETER :: dp   = KIND( 1.0D0 )
      REAL(kind=dp), parameter :: ZERO = 0.0_dp, ONE = 1.0_dp, half = 0.5_dp

!  subroutine arguments
!  --------------------

!  Dimensions

      INTEGER, INTENT(IN)         :: MAXLAYERS

!  Given layer index and Fourier number (inputs)

      INTEGER, INTENT(IN)         :: N
      INTEGER, INTENT(IN)         :: FOURIER

!  Stream value

      REAL(kind=dp), INTENT(IN)   :: STREAM_VALUE

!  Polynomials

      REAL(kind=dp), INTENT(IN)   :: PXSQ
      
!  OMEGA and ASYMM

      REAL(kind=dp), INTENT(IN)   :: OMEGA ( MAXLAYERS )
      REAL(kind=dp), INTENT(IN)   :: ASYMM ( MAXLAYERS )

!  optical thickness

      REAL(kind=dp), INTENT(IN)   :: DELTAU_VERT(MAXLAYERS)

!  Solutions to the homogeneous RT equations 
!  -----------------------------------------

!  local matrices for eigenvalue computation

      REAL(kind=dp), INTENT(INOUT)  :: SAB(MAXLAYERS), DAB(MAXLAYERS)

!  Eigensolutions

      REAL(kind=dp), INTENT(INOUT)  :: EIGENVALUE(MAXLAYERS)
      REAL(kind=dp), INTENT(INOUT)  :: EIGENTRANS(MAXLAYERS)

!  UP and down solutions

      REAL(kind=dp), INTENT(INOUT)  :: XPOS(2,MAXLAYERS)
      REAL(kind=dp), INTENT(INOUT)  :: XNEG(2,MAXLAYERS)

!  Gamma, Alpha and Epsilon

      REAL(kind=dp), INTENT(INOUT)  :: GAMMA(MAXLAYERS)
      REAL(kind=dp), INTENT(INOUT)  :: ALPHA(MAXLAYERS)
      REAL(kind=dp), INTENT(INOUT)  :: EPSON(MAXLAYERS)

!  Green;s function normalization factors
!    Introduced for [V2p3, Mark 10]

      REAL(kind=dp), INTENT(INOUT)  :: NORM_SAVED(MAXLAYERS)

!  Local variables
!  ---------------

!    parameter Introduced for [V2p3, Mark 10]

      REAL(kind=dp) :: EP, EM, XINV, OMEGA_ASYMM_3, DIFVEC, HELP,OG3P,y,b,yob
      REAL(kind=dp), parameter :: MAX_TAU_QPATH = 200.0_dp

!  Develop Sum and Difference matrices, set Eigenvalue, Alpha, Epsilon

      XINV = one / STREAM_VALUE
      OMEGA_ASYMM_3 = 3.0_dp * OMEGA(N) * ASYMM(N)
      OG3P = PXSQ * OMEGA_ASYMM_3
      ALPHA(N) = SQRT ( one - OG3P )
      EPSON(N) = SQRT ( one - OMEGA(N) )
      if ( fourier.eq.0) then
        EP = OMEGA(N) + OG3P
        EM = OMEGA(N) - OG3P
      Else if ( fourier .eq. 1 ) then
        EP = OG3P
        EM = OG3P
      ENDIF
      SAB(N) = XINV * ( ( EP + EM ) * 0.5_dp - one )
      DAB(N) = XINV * ( ( EP - EM ) * 0.5_dp - one )
      EIGENVALUE(N) = SQRT(SAB(N)*DAB(N))

!  Eigentrans, defined properly. [V2p3, Mark 10]

      HELP = EIGENVALUE(N)*DELTAU_VERT(N)
      IF ( HELP .GT. MAX_TAU_QPATH ) THEN
         EIGENTRANS(N) = zero
      ELSE
         EIGENTRANS(N) = EXP(-HELP)
      ENDIF

!write(*,*)n,eigenvalue(n),-dab(n),-sab(n) ! pause

!  Auxiliary equation to get up and down solutions

      DIFVEC = - SAB(N) / EIGENVALUE(N)
      XPOS(1,N) = 0.5d0 * ( one + DIFVEC )
      XPOS(2,N) = 0.5d0 * ( one - DIFVEC )

!  Symmetry and Gamma

      XNEG(1,N) = XPOS(2,N)
      XNEG(2,N) = XPOS(1,N)
      GAMMA(N)  = XPOS(2,N) / XPOS(1,N)

!  Taylor series on Gamma, conservative case

      if ( FOURIER.eq.0.and.EPSON(N).lt.1.0e-04 ) then
         y = EPSON(N) ; b = ALPHA(N) ; yob = y/b
         GAMMA(N) =  one-2.0_dp*yob*(one-yob*(one-yob*(one-yob*(one-yob*(one-yob)))))
      endif

!  Taylor series on Gamma, Pure absoprtion case

      if ( FOURIER.eq.0 .and. OMEGA(N).lt.1.0e-06 ) then
         b = 3.0_dp * pxsq * asymm(n)
         GAMMA(N) = 0.25_dp * omega(n) * (one-b) * ( one + half*omega(n)*(one+b) )
      endif

!  Green's function norm
!    Introduced for Version 2p4, Mark 10

      NORM_SAVED(N) = STREAM_VALUE *  &
              ( XPOS(1,N)*XPOS(1,N) - XPOS(2,N)*XPOS(2,N) )

!        write(*,*)n,-stream_value * SAB(n) / eigenvalue(n), NORM_SAVED(N) 

!  debug
!      if (fourier.eq.0.and.n.gt.50)write(45,'(i4,1p3e24.12)')n,EIGENTRANS(N),norm_saved(n),EIGENVALUE(N)
!      if ( n.eq.23)pause

!  Finish

      RETURN
END SUBROUTINE TWOSTREAM_HOM_SOLUTION

!

SUBROUTINE TWOSTREAM_HOM_USERSOLUTION_PD &
         ( MAXLAYERS, MAX_USER_STREAMS, CONSSCAT,          & ! Dimensions/CSFlag
           N_USER_STREAMS, N, FOURIER, STREAM_VALUE, PX11, & ! Input
           USER_STREAMS, ULP, XPOS, OMEGA, ASYMM,          & ! Input
           U_XPOS, U_XNEG, U_HELP_P, U_HELP_M )              ! Output

!  4/12/24. Add consscat flag

      implicit none

!  precision and parameters

      INTEGER      , PARAMETER :: dp   = KIND( 1.0D0 )
      REAL(kind=dp), parameter :: ZERO = 0.0_dp, ONE = 1.0_dp

!  subroutine input arguments
!  --------------------------

!  Dimensions

      INTEGER, INTENT(IN)        :: MAXLAYERS, MAX_USER_STREAMS

!  4/12/24. Flag for conservative scattering.

      LOGICAL, intent(in)        :: CONSSCAT

!  Numbers

      INTEGER, INTENT(IN)        :: N_USER_STREAMS

!  Given layer index and Fourier number (inputs)

      INTEGER, INTENT(IN)        :: N, FOURIER

!  Stream value and polynomial

      REAL(kind=dp), INTENT(IN)  :: STREAM_VALUE, PX11

!  User-defined post-processing stream directions

      REAL(kind=dp), INTENT(IN)  :: USER_STREAMS ( MAX_USER_STREAMS )
      REAL(kind=dp), INTENT(IN)  :: ULP          ( MAX_USER_STREAMS )

!  OMEGA and ASYMM

      REAL(kind=dp), INTENT(IN)  :: OMEGA ( MAXLAYERS )
      REAL(kind=dp), INTENT(IN)  :: ASYMM ( MAXLAYERS )

!  UP and down solutions

      REAL(kind=dp), INTENT(IN)  :: XPOS(2,MAXLAYERS)

!  Subroutine output arguments
!  ---------------------------

!  Eigenvectors defined at user-defined stream angles

      REAL(kind=dp), INTENT(INOUT) ::  U_XPOS(MAX_USER_STREAMS,MAXLAYERS)
      REAL(kind=dp), INTENT(INOUT) ::  U_XNEG(MAX_USER_STREAMS,MAXLAYERS)

      REAL(kind=dp), INTENT(OUT) ::    U_HELP_P(0:1)
      REAL(kind=dp), INTENT(OUT) ::    U_HELP_M(0:1)

!  Local variables
!  ---------------

      INTEGER       :: UM
      REAL(kind=dp) :: SUM_NEG, SUM_POS
      REAL(kind=dp) :: OMEGA_MOM, HMU_STREAM

!  zero the user solutions
!  4/12/24. return if conservative scattering

      DO UM = 1, N_USER_STREAMS
        U_XPOS(UM,N) = zero
        U_XNEG(UM,N) = zero
      ENDDO
      IF ( CONSSCAT ) RETURN

!  Eigenvector interpolation to user-defined angles
!  ------------------------------------------------

!  For each moment, do inner sum over computational angles
!  for the positive and negative eigenvectors

      HMU_STREAM = 0.5d0 * STREAM_VALUE
      if ( fourier.eq.0) then
        u_help_p(0) = ( XPOS(2,N) + XPOS(1,N) ) * 0.5d0
        u_help_p(1) = ( XPOS(2,N) - XPOS(1,N) ) * HMU_STREAM
        u_help_M(0) =   u_help_p(0)
        u_help_M(1) = - u_help_p(1)
      else
        u_help_p(1) = - ( XPOS(2,N) + XPOS(1,N) ) * PX11 * 0.5d0
        u_help_M(1) = u_help_p(1)
      endif

!  Now sum over harmonic contributions at each user-defined stream

      OMEGA_MOM = 3.0d0 * OMEGA(N) * ASYMM(N)
      DO UM = 1, N_USER_STREAMS
        if (fourier.eq.0 ) then
          sum_pos = u_help_p(0) * omega(N) &
                 +  u_help_p(1) * omega_mom * user_streams(um)
          sum_neg = u_help_m(0) * omega(N) &
                 +  u_help_m(1) * omega_mom * user_streams(um)
        else
          sum_pos = u_help_p(1) * omega_mom * ulp(um)
          sum_neg = u_help_m(1) * omega_mom * ulp(um)
        endif
        U_XPOS(UM,N) = SUM_POS
        U_XNEG(UM,N) = SUM_NEG
      ENDDO

!  debug
!      if (fourier.eq.1)
!     &   write(57,'(i4,1p2e24.12)')n,u_xpos(1,N),u_xneg(1,N)

!  Finish

      RETURN
END SUBROUTINE TWOSTREAM_HOM_USERSOLUTION_PD

!

SUBROUTINE TWOSTREAM_HOM_USERSOLUTION_TD &
         ( MAXLAYERS, MAX_USER_STREAMS, CONSSCAT,          & ! Dimensions/CSFlag
           N_USER_STREAMS, N, FOURIER, STREAM_VALUE, PX11, & ! Input
           USER_STREAMS, ULP, XPOS, OMEGA, ASYMM, GAMMA,   & ! Input
           U_XPOS, U_XNEG, U_HELP_P, U_HELP_M )              ! Output

!  4/12/24. Add consscat flag

      implicit none

!  precision and parameters

      INTEGER      , PARAMETER :: dp   = KIND( 1.0D0 )
      REAL(kind=dp), parameter :: ZERO = 0.0_dp, ONE = 1.0_dp

!  subroutine input arguments
!  --------------------------

!  Dimensions

      INTEGER, INTENT(IN)        :: MAXLAYERS, MAX_USER_STREAMS

!  4/12/24. Flag for conservative scattering.

      LOGICAL, intent(in)        :: CONSSCAT

!  Numbers

      INTEGER, INTENT(IN)        :: N_USER_STREAMS

!  Given layer index and Fourier number (inputs)

      INTEGER, INTENT(IN)        :: N, FOURIER

!  Stream value and polynomial

      REAL(kind=dp), INTENT(IN)  :: STREAM_VALUE, PX11

!  User-defined post-processing stream directions

      REAL(kind=dp), INTENT(IN)  :: USER_STREAMS ( MAX_USER_STREAMS )
      REAL(kind=dp), INTENT(IN)  :: ULP          ( MAX_USER_STREAMS )

!  OMEGA and ASYMM, and GAMMA

      REAL(kind=dp), INTENT(IN)  :: OMEGA ( MAXLAYERS )
      REAL(kind=dp), INTENT(IN)  :: ASYMM ( MAXLAYERS )
      REAL(kind=dp), INTENT(IN)  :: GAMMA ( MAXLAYERS )

!  UP and down solutions

      REAL(kind=dp), INTENT(IN)  :: XPOS(2,MAXLAYERS)

!  Subroutine output arguments
!  ---------------------------

!  Eigenvectors defined at user-defined stream angles

      REAL(kind=dp), INTENT(INOUT) ::  U_XPOS(MAX_USER_STREAMS,MAXLAYERS)
      REAL(kind=dp), INTENT(INOUT) ::  U_XNEG(MAX_USER_STREAMS,MAXLAYERS)

      REAL(kind=dp), INTENT(OUT) ::    U_HELP_P(0:1)
      REAL(kind=dp), INTENT(OUT) ::    U_HELP_M(0:1)

!  Local variables
!  ---------------

      INTEGER       :: UM
      REAL(kind=dp) :: SUM_NEG, SUM_POS
      REAL(kind=dp) :: OMEGA_MOM, HMU_STREAM

!  zero the user solutions
!  4/12/24. return if conservative scattering

      DO UM = 1, N_USER_STREAMS
        U_XPOS(UM,N) = zero
        U_XNEG(UM,N) = zero
      ENDDO
      IF ( CONSSCAT ) RETURN

!  Eigenvector interpolation to user-defined angles
!  ------------------------------------------------

!  For each moment, do inner sum over computational angles
!  for the positive and negative eigenvectors

      HMU_STREAM = 0.5d0 * STREAM_VALUE
      if ( fourier.eq.0) then
!        u_help_p(0) = 0.5d0 * ( one + GAMMA(N) )
!        u_help_p(1) = - HMU_STREAM * ( one - GAMMA(N) )
        u_help_p(0) = 0.5d0 
        u_help_p(1) = HMU_STREAM 
        u_help_M(0) =   u_help_p(0)
        u_help_M(1) = - u_help_p(1)
      else
        u_help_p(1) = - PX11 * 0.5d0 * ( one + GAMMA(N) )
        u_help_M(1) = u_help_p(1)
      endif

!  Now sum over harmonic contributions at each user-defined stream

      OMEGA_MOM = 3.0d0 * OMEGA(N) * ASYMM(N)
      DO UM = 1, N_USER_STREAMS
        if (fourier.eq.0 ) then
          sum_pos = u_help_p(0) * omega(N) &
                 +  u_help_p(1) * omega_mom * user_streams(um)
          sum_neg = u_help_m(0) * omega(N) &
                 +  u_help_m(1) * omega_mom * user_streams(um)
        else
          sum_pos = u_help_p(1) * omega_mom * ulp(um)
          sum_neg = u_help_m(1) * omega_mom * ulp(um)
        endif
        U_XPOS(UM,N) = SUM_POS
        U_XNEG(UM,N) = SUM_NEG
      ENDDO

!  debug
!      if (fourier.eq.1)
!     &   write(57,'(i4,1p2e24.12)')n,u_xpos(1,N),u_xneg(1,N)

!  Finish

      RETURN
END SUBROUTINE TWOSTREAM_HOM_USERSOLUTION_TD

!

SUBROUTINE TWOSTREAM_HMULT_MASTER &
           ( MAXLAYERS, MAX_USER_STREAMS, CONSSCAT,         & ! Dimensions
             TAYLOR_ORDER, TAYLOR_SMALL, DELTAUS, ASYMM,    & ! Inputs 
             NLAYERS, N_USER_STREAMS, USER_SECANTS, STREAM_VALUE,  & ! Inputs
             EIGENVALUE, EIGENTRANS, T_DELT_USERM, TGMFUNC,        & ! Inputs
             ZETA_M, ZETA_P, HMULT_1, HMULT_2, CSPPHOMUP, CSPPHOMDN )             ! Output

      implicit none

!  precision and parameters

      INTEGER      , PARAMETER :: dp   = KIND( 1.0D0 )
      REAL(kind=dp), parameter :: ZERO = 0.0_dp, ONE = 1.0_dp

!  Input arguments
!  ===============

!  Dimensions

      INTEGER, INTENT(IN)        :: MAXLAYERS, MAX_USER_STREAMS

!  Order of Taylor series (including terms up to EPS^n).
!    Introduced for [V2p3, Mark 10]

      INTEGER      , intent(in)  :: TAYLOR_ORDER
      REAL(kind=dp), INTENT(IN)  :: TAYLOR_SMALL

!  Input Optical depths required for Taylor-series limiting cases
!    Introduced for [V2p3, Mark 10]

      REAL(kind=dp), intent(in)  :: DELTAUS(MAXLAYERS)
      REAL(kind=dp), intent(in)  :: ASYMM  (MAXLAYERS)

!  Numbers

      INTEGER, INTENT(IN)        :: NLAYERS, N_USER_STREAMS

!  User secants (formerly streams). [V2p3, Mark 10]

      REAL(kind=dp), INTENT(IN)  :: USER_SECANTS ( MAX_USER_STREAMS )

!  Stream Value introduced, 3/28/18

      REAL(kind=dp), INTENT(IN)  :: STREAM_VALUE

!  Transmittance factors for user-defined stream angles

      REAL(kind=dp), INTENT(IN)  :: T_DELT_USERM ( MAXLAYERS, MAX_USER_STREAMS )

!  Eigensolutions

      REAL(kind=dp), INTENT(IN)  :: EIGENVALUE(MAXLAYERS)
      REAL(kind=dp), INTENT(IN)  :: EIGENTRANS(MAXLAYERS)

!  Conservative scattering solutions

      LOGICAL      , INTENT(IN)  :: CONSSCAT(MAXLAYERS)
      REAL(kind=dp), INTENT(IN)  :: TGMFUNC(MAXLAYERS)

!  Output = Global multipliers
!  ===========================

!  coefficient functions for user-defined angles.
!    Formerly, defined as inverses.  [V2p3, Mark 10]

      REAL(kind=dp), INTENT(INOUT) :: ZETA_P(MAX_USER_STREAMS,MAXLAYERS)
      REAL(kind=dp), INTENT(INOUT) :: ZETA_M(MAX_USER_STREAMS,MAXLAYERS)

!  Integrated homogeneous solution multipliers, whole layer

      REAL(kind=dp), INTENT(INOUT) :: HMULT_1(MAX_USER_STREAMS,MAXLAYERS)
      REAL(kind=dp), INTENT(INOUT) :: HMULT_2(MAX_USER_STREAMS,MAXLAYERS)

!  outputs (3/28/18 programmed)

      REAL(kind=dp), INTENT(INOUT) :: CSPPHOMUP(2,MAX_USER_STREAMS,MAXLAYERS)
      REAL(kind=dp), INTENT(INOUT) :: CSPPHOMDN(2,MAX_USER_STREAMS,MAXLAYERS)

!  Local variables
!  ---------------

      INTEGER       :: UM, N
      REAL(kind=dp) :: UDEL, SM, EPS, ZDEL, ZUDEL
      REAL(kind=dp) :: K0, K1, UTGP, UTGM, SMINV

!  whole layer multipliers
!  -----------------------

!  Start loops over layers and user-streams
!    Only done if layers are flagged

      DO N = 1, NLAYERS
        IF ( CONSSCAT(N) ) THEN
          DO UM = 1, N_USER_STREAMS
            UDEL = T_DELT_USERM(N,UM)
            SM   = USER_SECANTS(UM) ; SMINV = ONE / SM
            K0 =   SMINV * ( ONE - UDEL )
            K1 = - SMINV * ( DELTAUS(N) * UDEL - K0 )
            UTGP = 0.5d0 * ( one + 3.0_dp * ASYMM(N) * SMINV * STREAM_VALUE )
            UTGM = 0.5d0 * ( one - 3.0_dp * ASYMM(N) * SMINV * STREAM_VALUE )
            CSPPHOMUP(1,UM,N) = SM * ( UTGM * K0 - TGMFUNC(N) * K1 )
            CSPPHOMUP(2,UM,N) = SM * ( UTGP * K0 + TGMFUNC(N) * K1 )
            CSPPHOMDN(1,UM,N) = SM * ( UTGP * K0 - TGMFUNC(N) * K1 )
            CSPPHOMDN(2,UM,N) = SM * ( UTGM * K0 + TGMFUNC(N) * K1 )
          ENDDO
        ELSE
          DO UM = 1, N_USER_STREAMS
            UDEL = T_DELT_USERM(N,UM)
            SM   = USER_SECANTS(UM)
            ZETA_P(UM,N) = SM + EIGENVALUE(N)
            ZETA_M(UM,N) = SM - EIGENVALUE(N)
            ZDEL    = EIGENTRANS(N)
            ZUDEL   = ZDEL * UDEL
            HMULT_2(UM,N) = SM * ( ONE - ZUDEL ) / ZETA_P(UM,N)
            IF ( ABS(ZETA_M(UM,N)) .LT. TAYLOR_SMALL ) THEN
              EPS = ZETA_M(UM,N)
              CALL TWOSTREAM_TAYLOR_SERIES_1 ( TAYLOR_ORDER, EPS, DELTAUS(N), UDEL, SM, HMULT_1(UM,N) )
            ELSE
              HMULT_1(UM,N) = SM * ( ZDEL - UDEL ) / ZETA_M(UM,N)
            ENDIF
          ENDDO
        ENDIF
      ENDDO

!  debug
!      do n = 1, 3
!        write(*,*)HMULT_1(1,N),HMULT_2(2,N)
!      enddo
!      pause

!  Finish

      RETURN
END SUBROUTINE TWOSTREAM_HMULT_MASTER

!

SUBROUTINE TWOSTREAM_GBEAM_SOLUTION &
         ( MAXLAYERS, MAXBEAMS,                                & ! Dimensions
           TAYLOR_ORDER, TAYLOR_SMALL, DELTAUS,                & ! Inputs 
           N, FOURIER, IBEAM, PI4, FLUX_FACTOR,                & ! Inputs
           BEAM_CUTOFF, PX0X, OMEGA, ASYMM,               & ! Inputs
           AVERAGE_SECANT, INITIAL_TRANS, T_DELT_MUBAR,        & ! Inputs
           XPOS, EIGENVALUE, EIGENTRANS, NORM_SAVED,           & ! Inputs
           GAMMA_M, GAMMA_P, DMI, DPI, ATERM_SAVE, BTERM_SAVE, & ! Output
           CFUNC, DFUNC, GFUNC_UP, GFUNC_DN, WUPPER, WLOWER )    ! Output

      implicit none

!  precision and parameters

      INTEGER      , PARAMETER :: dp   = KIND( 1.0D0 )
      REAL(kind=dp), parameter :: ZERO = 0.0_dp, ONE = 1.0_dp

!  subroutine arguments
!  --------------------

!  Dimensions

      INTEGER, INTENT(IN)         :: MAXLAYERS, MAXBEAMS

!  Order of Taylor series (including terms up to EPS^n).
!    Introduced for [V2p3, Mark 10]

      INTEGER      , intent(in)  :: TAYLOR_ORDER
      REAL(kind=dp), INTENT(IN)  :: TAYLOR_SMALL

!  Input Optical depths required for Taylor-series limiting cases
!    Introduced for [V2p3, Mark 10]

      REAL(kind=dp), intent(in)  :: DELTAUS(MAXLAYERS)

!  Given layer index and Fourier number, Beam number (inputs)

      INTEGER, INTENT(IN)         :: N, FOURIER, IBEAM

!  Flux factor

      REAL(kind=dp), INTENT(IN)   :: FLUX_FACTOR, PI4

!  Last layer to include Particular integral solution

      INTEGER, INTENT(IN)         :: BEAM_CUTOFF(MAXBEAMS)

!  Beam SZA polynomial factors

      REAL(kind=dp), INTENT(IN)   :: PX0X(MAXBEAMS)

!  OMEGA and ASYMM

      REAL(kind=dp), INTENT(IN)   :: OMEGA ( MAXLAYERS )
      REAL(kind=dp), INTENT(IN)   :: ASYMM ( MAXLAYERS )

!  Average-secant and initial tramsittance factors for solar beams.

      REAL(kind=dp), INTENT(IN)   :: INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS )
      REAL(kind=dp), INTENT(IN)   :: AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS )
      REAL(kind=dp), INTENT(IN)   :: T_DELT_MUBAR   ( MAXLAYERS, MAXBEAMS )

!  UP and down solutions

      REAL(kind=dp), INTENT(IN)  :: XPOS(2,MAXLAYERS)

!  Eigenvalues and eigentransmittance

      REAL(kind=dp), INTENT(IN)   :: EIGENVALUE(MAXLAYERS)
      REAL(kind=dp), INTENT(IN)   :: EIGENTRANS(MAXLAYERS)

!  Green;s function normalization factors
!    Introduced for [V2p3, Mark 10]

      REAL(kind=dp), INTENT(IN)  :: NORM_SAVED(MAXLAYERS)

!  subroutine output arguments
!  ===========================

!mick fix 6/29/11 - change most outputs from "out" to "inout"

!  Saved quantities for the Green function solution

      REAL(kind=dp), intent(inout) :: ATERM_SAVE(MAXLAYERS)
      REAL(kind=dp), intent(inout) :: BTERM_SAVE(MAXLAYERS)
      REAL(kind=dp), intent(inout) :: DMI, DPI

!  Layer C and D functions

      REAL(kind=dp), intent(inout) :: CFUNC(MAXLAYERS)
      REAL(kind=dp), intent(inout) :: DFUNC(MAXLAYERS)

!  Green function Multipliers for solution

      REAL(kind=dp), intent(inout) :: GFUNC_UP(MAXLAYERS)
      REAL(kind=dp), intent(inout) :: GFUNC_DN(MAXLAYERS)

!  Holding arrays for Multiplier coefficients

      REAL(kind=dp), intent(inout) :: GAMMA_M(MAXLAYERS)
      REAL(kind=dp), intent(inout) :: GAMMA_P(MAXLAYERS)

!  Solutions at layer boundaries

      REAL(kind=dp), INTENT(INOUT)  :: WUPPER(2,MAXLAYERS)
      REAL(kind=dp), INTENT(INOUT)  :: WLOWER(2,MAXLAYERS)

!  Help variables
!  --------------

      INTEGER       :: I
      REAL(kind=dp) :: TP, TM, EPS, CONST, SECBAR
      REAL(kind=dp) :: WDEL, ZDEL, ZWDEL, F1, SUM_LA, SUM_LB, OMEGA_ASYMM

!  Flux factor

      F1 = FLUX_FACTOR / PI4

!  No particular solution beyond the cutoff layer
!  Or no scattering in this layer...
!  ... Zero the boundary layer values and exit

      IF ( N .GT. BEAM_CUTOFF(IBEAM) ) THEN
        DO I = 1, 2
          WUPPER(I,N) = zero
          WLOWER(I,N) = zero
        ENDDO
        RETURN
      ENDIF

!  constants for the layer

      SECBAR = AVERAGE_SECANT(N,IBEAM)
      CONST  = INITIAL_TRANS (N,IBEAM)
      WDEL   = T_DELT_MUBAR  (N,IBEAM)

!  Optical depth integrations for the discrete ordinate solution
!  =============================================================

      GAMMA_P(N) = SECBAR + EIGENVALUE(N)
      GAMMA_M(N) = SECBAR - EIGENVALUE(N)
      ZDEL  = EIGENTRANS(N)
      ZWDEL = ZDEL * WDEL
      IF ( ABS(GAMMA_M(N)) .LT. TAYLOR_SMALL ) THEN
         EPS = GAMMA_M(N)
         CALL TWOSTREAM_TAYLOR_SERIES_1 ( TAYLOR_ORDER, EPS, DELTAUS(N), WDEL, ONE, CFUNC(N) )
      ELSE
         CFUNC(N) =  ( ZDEL - WDEL ) / GAMMA_M(N)
      ENDIF
      DFUNC(N)  = ( one - ZWDEL ) / GAMMA_P(N)

!  Help quantitiesfor Green's function

      OMEGA_ASYMM = OMEGA(N) * ASYMM(N) * 3.0_dp
      if ( fourier.eq.0) then
         TP = OMEGA(N) + PX0X(IBEAM) * OMEGA_ASYMM
         TM = OMEGA(N) - PX0X(IBEAM) * OMEGA_ASYMM
      Else if ( fourier .eq. 1 ) then
         TP = PX0X(IBEAM) * OMEGA_ASYMM
         TM = PX0X(IBEAM) * OMEGA_ASYMM
      ENDIF
      DPI = TP * F1
      DMI = TM * F1

!  Green function multipliers GFUNC

      SUM_LA  = DPI*XPOS(1,N)+DMI*XPOS(2,N)
      SUM_LB  = DMI*XPOS(1,N)+DPI*XPOS(2,N)
      ATERM_SAVE(N) = SUM_LA / NORM_SAVED(N)
      BTERM_SAVE(N) = SUM_LB / NORM_SAVED(N)
      GFUNC_DN(N) = CFUNC(N) * ATERM_SAVE(N) * CONST
      GFUNC_UP(N) = DFUNC(N) * BTERM_SAVE(N) * CONST

!  particular integrals at lower and upper boundaries

      WUPPER(1,N) = GFUNC_UP(N)*XPOS(2,N)
      WUPPER(2,N) = GFUNC_UP(N)*XPOS(1,N)
      WLOWER(1,N) = GFUNC_DN(N)*XPOS(1,N)
      WLOWER(2,N) = GFUNC_DN(N)*XPOS(2,N)

!      if ( FOURIER.eq.0)write(*,*)'Greens',Fourier, ibeam, n,WUPPER(1:2,n)

!  Finish

      RETURN
END SUBROUTINE TWOSTREAM_GBEAM_SOLUTION

!

SUBROUTINE TWOSTREAM_CBEAM_SOLUTION &
         ( MAXLAYERS, MAXBEAMS, N, FOURIER, IBEAM, PI4,  & ! Inputs
           FLUX_FACTOR, BEAM_CUTOFF, STREAM_VALUE, PX0X, & ! Inputs
           AVERAGE_SECANT, INITIAL_TRANS, T_DELT_MUBAR,  & ! Inputs
           OMEGA, ASYMM, SAB, DAB, EIGENVALUE,           & ! Inputs
           QSUMVEC, QDIFVEC, QVEC, WVEC, WUPPER, WLOWER )  ! In/Out

      implicit none

!  precision and parameters

      INTEGER      , PARAMETER :: dp   = KIND( 1.0D0 )
      REAL(kind=dp), parameter :: ZERO = 0.0_dp, ONE = 1.0_dp

!  subroutine arguments
!  --------------------

!  Dimensions

      INTEGER, INTENT(IN)         :: MAXLAYERS, MAXBEAMS

!  Given layer index and Fourier number, Beam number (inputs)

      INTEGER, INTENT(IN)         :: N, FOURIER, IBEAM

!  Flux factor

      REAL(kind=dp), INTENT(IN)   :: FLUX_FACTOR, PI4

!  Last layer to include Particular integral solution

      INTEGER, INTENT(IN)         :: BEAM_CUTOFF(MAXBEAMS)

!  Stream value

      REAL(kind=dp), INTENT(IN)   :: STREAM_VALUE

!  Beam SZA polynomial factors

      REAL(kind=dp), INTENT(IN)   :: PX0X(MAXBEAMS)

!  Average-secant and initial tramsittance factors for solar beams.

      REAL(kind=dp), INTENT(IN)   :: INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS )
      REAL(kind=dp), INTENT(IN)   :: AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS )
      REAL(kind=dp), INTENT(IN)   :: T_DELT_MUBAR   ( MAXLAYERS, MAXBEAMS )

!  OMEGA and ASYMM

      REAL(kind=dp), INTENT(IN)   :: OMEGA ( MAXLAYERS )
      REAL(kind=dp), INTENT(IN)   :: ASYMM ( MAXLAYERS )

!  local matrices from eigenvalue computation

      REAL(kind=dp), INTENT(IN)   :: SAB(MAXLAYERS), DAB(MAXLAYERS)

!  Eigenvalues

      REAL(kind=dp), INTENT(IN)   :: EIGENVALUE(MAXLAYERS)

!  Output variables
!  ----------------

!  Auxiliary vectors

      REAL(kind=dp), INTENT(INOUT)  :: QDIFVEC(MAXLAYERS)
      REAL(kind=dp), INTENT(INOUT)  :: QSUMVEC(MAXLAYERS)
      REAL(kind=dp), INTENT(INOUT)  :: QVEC   (MAXLAYERS)

!  Beam solution

      REAL(kind=dp), INTENT(INOUT)  :: WVEC(2,MAXLAYERS)

!  Solutions at layer boundaries

      REAL(kind=dp), INTENT(INOUT)  :: WUPPER(2,MAXLAYERS)
      REAL(kind=dp), INTENT(INOUT)  :: WLOWER(2,MAXLAYERS)

!  Help variables
!  --------------

      INTEGER       :: I
      REAL(kind=dp) :: TP, TM, INV_X0SQ, SECBAR, XINV, F1
      REAL(kind=dp) :: HELP, TRANS1, TRANS2
      REAL(kind=dp) :: QMAT, QDIF, OMEGA_ASYMM

!  Flux factor

      F1 = FLUX_FACTOR / PI4

!  No particular solution beyond the cutoff layer
!  ... Zero the boundary layer values and exit

      IF ( N .GT. BEAM_CUTOFF(IBEAM) ) THEN
        WUPPER(1:2,N) = zero
        WLOWER(1:2,N) = zero
        RETURN
      ENDIF

!  set local values

      SECBAR   = AVERAGE_SECANT(N,IBEAM)
      INV_X0SQ = SECBAR * SECBAR

!  Set up sum and differences for Beam source terms
!  ( sum may be required again in linearization )

      XINV = one / STREAM_VALUE
      OMEGA_ASYMM = OMEGA(N) * ASYMM(N) * 3.0_dp
      if ( fourier.eq.0) then
        TP = OMEGA(N) + PX0X(IBEAM) * OMEGA_ASYMM
        TM = OMEGA(N) - PX0X(IBEAM) * OMEGA_ASYMM
      Else if ( fourier .eq. 1 ) then
        TP = PX0X(IBEAM) * OMEGA_ASYMM
        TM = PX0X(IBEAM) * OMEGA_ASYMM
      ENDIF
      QSUMVEC(N) =  F1 * ( TP + TM ) * XINV
      QDIFVEC(N) =  F1 * ( TP - TM ) * XINV

!  The reduced problem: QMAT. W = QVEC (Overwrite QVEC)
!  Restore up and down solutions

      QMAT = EIGENVALUE(N) * EIGENVALUE(N) - INV_X0SQ
      HELP = - DAB(N) * QSUMVEC(N)
      QVEC(N) = HELP + QDIFVEC(N) * SECBAR
      QVEC(N) = QVEC(N) / QMAT
      HELP = - SAB(N) * QVEC(N)
      QDIF = ( HELP - QSUMVEC(N) ) / SECBAR
      WVEC(1,N) = 0.5_dp * ( QVEC(N) + QDIF )
      WVEC(2,N) = 0.5_dp * ( QVEC(N) - QDIF )

!  Values at the layer boundaries
!  (transmittance factors have been determined in SETUPS module)

      TRANS1 = INITIAL_TRANS(N,IBEAM)
      TRANS2 = T_DELT_MUBAR(N,IBEAM) * TRANS1

      DO I = 1, 2
        WUPPER(I,N) = WVEC(I,N)*TRANS1
        WLOWER(I,N) = WVEC(I,N)*TRANS2
      ENDDO

!  Finish

      RETURN
END SUBROUTINE TWOSTREAM_CBEAM_SOLUTION

!

SUBROUTINE TWOSTREAM_CBEAM_USERSOLUTION &
       ( DO_UPWELLING, DO_DNWELLING,                             & ! Input
         NLAYERS, NBEAMS, N_USER_STREAMS, N, FOURIER, IBEAM,     & ! Input
         FLUX_FACTOR, BEAM_CUTOFF, STREAM_VALUE, PI4, PX11, & ! Input
         OMEGA, ASYMM, USER_STREAMS, ULP, WVEC,                  & ! Input
         U_WPOS2, U_WNEG2, W_HELP )                                ! Output

      implicit none

!  precision

      INTEGER, PARAMETER :: dp     = KIND( 1.0D0 )
      REAL(kind=dp), parameter :: ZERO = 0.0_dp, ONE = 1.0_dp

!  subroutine arguments
!  --------------------

!  Flags

      LOGICAL, INTENT(IN)         :: DO_UPWELLING, DO_DNWELLING

!  Numbers

      INTEGER, INTENT(IN)         :: NLAYERS, NBEAMS, N_USER_STREAMS

!  Given layer index and Fourier number, Beam number (inputs)

      INTEGER, INTENT(IN)         :: N, FOURIER, IBEAM

!  Flux factor

      REAL(kind=dp), INTENT(IN)   :: FLUX_FACTOR, PI4

!  Last layer to include Particular integral solution

      INTEGER, INTENT(IN)         :: BEAM_CUTOFF(NBEAMS)

!  Stream value and polynomial

      REAL(kind=dp), INTENT(IN)   :: STREAM_VALUE, PX11

!  Beam SZA cosines. Not required for MS-mode only
!      REAL(kind=dp), INTENT(IN)   :: X0(NBEAMS)
!      REAL(kind=dp), INTENT(IN)   :: POX(NBEAMS)

!  OMEGA and ASYMM

      REAL(kind=dp), INTENT(IN)   :: OMEGA ( NLAYERS )
      REAL(kind=dp), INTENT(IN)   :: ASYMM ( NLAYERS )

!  User streams

      REAL(kind=dp), INTENT(IN)   :: USER_STREAMS ( N_USER_STREAMS )
      REAL(kind=dp), INTENT(IN)   :: ULP          ( N_USER_STREAMS )

!  Beam solution

      REAL(kind=dp), INTENT(IN)   :: WVEC(2,NLAYERS)

!  Subroutine output arguments
!  ---------------------------

!  Single-scatter Particular beam solutions at user-defined angles
!    NOT REQUIRED, MS_MODE only
!      REAL(kind=dp), INTENT(INOUT) :: U_WPOS1(N_USER_STREAMS,NLAYERS)
!      REAL(kind=dp), INTENT(INOUT) :: U_WNEG1(N_USER_STREAMS,NLAYERS)

!  Diffuse-term Particular beam solutions at user-defined angles

      REAL(kind=dp), INTENT(INOUT) :: U_WPOS2(N_USER_STREAMS,NLAYERS)
      REAL(kind=dp), INTENT(INOUT) :: U_WNEG2(N_USER_STREAMS,NLAYERS)

!  Saved help variables

      REAL(kind=dp), INTENT(OUT) ::   W_HELP(0:1)

!  Local variables
!  ---------------

      INTEGER       :: UM
      REAL(kind=dp) :: POS2, F1, OMEGA_MOM
      REAL(kind=dp) :: HELP2(0:1), HMU_STREAM

!  No particular solution beyond the cutoff layer
!  ... Zero the user solutions and exit

      IF ( N .GT. BEAM_CUTOFF(IBEAM) ) THEN
        U_WPOS2(1:N_USER_STREAMS,N) = zero
        U_WNEG2(1:N_USER_STREAMS,N) = zero
        RETURN
      ENDIF

!  Scattering solutions
!  ====================

!  Starter quantities

      F1 = FLUX_FACTOR / PI4
      OMEGA_MOM = 3.0d0 * OMEGA(N) * ASYMM(N)
      HMU_STREAM = STREAM_VALUE * 0.5d0

!  For each moment do inner sum over computational angles

      if ( fourier.eq.0) then
        w_help(0) = ( WVEC(2,N) + WVEC(1,N) ) * 0.5d0
        w_help(1) = ( WVEC(2,N) - WVEC(1,N) ) * HMU_STREAM
        help2(0)  = w_help(0) * omega(n)
        help2(1)  = w_help(1) * omega_mom
      else
        w_help(1) = - ( WVEC(2,N) + WVEC(1,N) ) * PX11 * 0.5d0
        help2(1)  = w_help(1) * omega_mom
      endif

!  Now sum over all harmonic contributions at each user-defined stream
!  Distinguish between upwelling and downwelling

      IF ( DO_UPWELLING ) THEN
        DO UM = 1, N_USER_STREAMS
          if (fourier.eq.0 ) then
            pos2 = help2(0) + help2(1)* user_streams(um)
          else
            pos2 = help2(1)* ULP(UM)
          endif
          U_WPOS2(UM,N) = POS2
        ENDDO
      ENDIF

      IF ( DO_DNWELLING ) THEN
        DO UM = 1, N_USER_STREAMS
          if (fourier.eq.0 ) then
            pos2 = help2(0) - help2(1)* user_streams(um)
          else
            pos2 = help2(1) * ulp(UM)
          endif
          U_WNEG2(UM,N) = POS2
        ENDDO
      ENDIF

! Finish

      RETURN
END SUBROUTINE TWOSTREAM_CBEAM_USERSOLUTION

!

SUBROUTINE TWOSTREAM_CONSSCAT_SOLUTION &
         ( MAXLAYERS, MAXBEAMS, DO_PLANPAR, N, IBEAM, PI4, MU0, & ! Inputs
           FLUX_FACTOR, BEAM_CUTOFF, STREAM_VALUE, ASYMM,  & ! Inputs
           INITIAL_TRANS, AVERAGE_SECANT, T_DELT_MUBAR,         & ! Inputs
           WVEC, WUPPER, WLOWER )                                 ! In/Out

      implicit none

!  precision and parameters

      INTEGER      , PARAMETER :: dp   = KIND( 1.0D0 )
      REAL(kind=dp), parameter :: ZERO = 0.0_dp, ONE = 1.0_dp

!  subroutine arguments
!  --------------------

!  Dimensions

      INTEGER, INTENT(IN)         :: MAXLAYERS, MAXBEAMS

!  Plane-parallel flag

      LOGICAL, INTENT(IN)         :: DO_PLANPAR

!  Given layer index and Beam number (inputs)

      INTEGER, INTENT(IN)         :: N, IBEAM

!  Solar angle cosine

      REAL(kind=dp), INTENT(IN)   :: MU0

!  Flux factor

      REAL(kind=dp), INTENT(IN)   :: FLUX_FACTOR, PI4

!  Last layer to include Particular integral solution

      INTEGER, INTENT(IN)         :: BEAM_CUTOFF(MAXBEAMS)

!  Stream value

      REAL(kind=dp), INTENT(IN)   :: STREAM_VALUE

!  ASYMM

      REAL(kind=dp), INTENT(IN)   :: ASYMM ( MAXLAYERS )

!  Average-secant and tramsittance factors for solar beams.

      REAL(kind=dp), INTENT(IN)   :: AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS )
      REAL(kind=dp), INTENT(IN)   :: INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS )
      REAL(kind=dp), INTENT(IN)   :: T_DELT_MUBAR   ( MAXLAYERS, MAXBEAMS )

!  Output variables
!  ----------------

!  Beam solution

      REAL(kind=dp), INTENT(INOUT)  :: WVEC(2,MAXLAYERS)

!  Solutions at layer boundaries

      REAL(kind=dp), INTENT(INOUT)  :: WUPPER(2,MAXLAYERS)
      REAL(kind=dp), INTENT(INOUT)  :: WLOWER(2,MAXLAYERS)

!  Help variables
!  --------------

      INTEGER       :: I
      REAL(kind=dp) :: F1, TRANS1, TRANS2, HELP, XINV, SECBAR, INV_X0SQ, TERM2

!  Flux factor

      F1 = FLUX_FACTOR / PI4

!  No particular solution beyond the cutoff layer
!  ... Zero the boundary layer values and exit

      IF ( N .GT. BEAM_CUTOFF(IBEAM) ) THEN
        WUPPER(1:2,N) = zero
        WLOWER(1:2,N) = zero
        RETURN
      ENDIF

!  set WVEC.
!    pseudo-spherical implemented, 3/28/18

      IF ( DO_PLANPAR ) THEN
        XINV = one / STREAM_VALUE
        HELP = MU0 * XINV
        WVEC(1,N) = - F1 * HELP * ( ONE + HELP ) 
        WVEC(2,N) =   F1 * HELP * ( ONE - HELP ) 
      ELSE
        XINV = one / STREAM_VALUE
        SECBAR   = AVERAGE_SECANT(N,IBEAM)
        INV_X0SQ = SECBAR * SECBAR
        HELP = XINV / SECBAR
        TERM2 = 3.0_dp * ASYMM(N) * ( SECBAR * MU0 - one ) / INV_X0SQ
        WVEC(1,N) = - F1 * ( HELP * ( ONE + HELP ) + TERM2 ) 
        WVEC(2,N) =   F1 * ( HELP * ( ONE - HELP ) - TERM2 )
      ENDIF

!  Values at the layer boundaries

      TRANS1 = INITIAL_TRANS(N,IBEAM)
      TRANS2 = T_DELT_MUBAR(N,IBEAM) * TRANS1
      DO I = 1, 2
        WUPPER(I,N) = WVEC(I,N)*TRANS1
        WLOWER(I,N) = WVEC(I,N)*TRANS2
      ENDDO

!  Finish

      RETURN
END SUBROUTINE TWOSTREAM_CONSSCAT_SOLUTION

!

SUBROUTINE EDDINGTON_BEAM_SOLUTION &
         ( MAXLAYERS, MAXBEAMS, N, IBEAM, PI4, MU0, FLUX_FACTOR,     & ! Inputs
           BEAM_CUTOFF, AVERAGE_SECANT, INITIAL_TRANS, T_DELT_MUBAR, & ! Inputs
           OMEGA_TOTAL, ASYMM_TOTAL, EIGENVALUE, GAMMA1, GAMMA2,     & ! Inputs
           WVEC, WUPPER, WLOWER )                                      ! In/Out

      implicit none

!  precision and parameters

      INTEGER      , PARAMETER :: dp   = KIND( 1.0D0 )
      REAL(kind=dp), parameter :: ZERO = 0.0_dp, ONE = 1.0_dp

!  subroutine arguments
!  --------------------

!  Dimensions

      INTEGER, INTENT(IN)         :: MAXLAYERS, MAXBEAMS

!  Given layer index and Fourier number, Beam number (inputs)

      INTEGER, INTENT(IN)         :: N, IBEAM

!  Flux factor

      REAL(kind=dp), INTENT(IN)   :: FLUX_FACTOR, PI4

!  SZA values

      REAL(kind=dp), INTENT(IN)   :: MU0 ( MAXBEAMS )

!  Last layer to include Particular integral solution

      INTEGER, INTENT(IN)         :: BEAM_CUTOFF(MAXBEAMS)

!  Average-secant and initial tramsittance factors for solar beams.

      REAL(kind=dp), INTENT(IN)   :: INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS )
      REAL(kind=dp), INTENT(IN)   :: AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS )
      REAL(kind=dp), INTENT(IN)   :: T_DELT_MUBAR   ( MAXLAYERS, MAXBEAMS )

!  OMEGA and ASYMM

      REAL(kind=dp), INTENT(IN)   :: OMEGA_TOTAL ( MAXLAYERS )
      REAL(kind=dp), INTENT(IN)   :: ASYMM_TOTAL ( MAXLAYERS )

!  Gamma1, Gamma2, Eigenvalues

      REAL(kind=dp), INTENT(IN)   :: GAMMA1(MAXLAYERS)
      REAL(kind=dp), INTENT(IN)   :: GAMMA2(MAXLAYERS)
      REAL(kind=dp), INTENT(IN)   :: EIGENVALUE(MAXLAYERS)

!  Output variables
!  ----------------

!  Beam solution

      REAL(kind=dp), INTENT(INOUT)  :: WVEC(2,MAXLAYERS)

!  Solutions at layer boundaries

      REAL(kind=dp), INTENT(INOUT)  :: WUPPER(2,MAXLAYERS)
      REAL(kind=dp), INTENT(INOUT)  :: WLOWER(2,MAXLAYERS)

!  Help variables
!  --------------

      REAL(kind=dp) :: MU0_INV, INV_X0SQ, SECBAR, FACTOR, F1
      REAL(kind=dp) :: GAMMA3, GAMMA4, QMAT, TRANS1, TRANS2

!  Flux factor. 6/29/24. Don't divide by PI4

!      F1 = FLUX_FACTOR / PI4
      F1 = FLUX_FACTOR

!  No particular solution beyond the cutoff layer
!  ... Zero the boundary layer values and exit

      IF ( N .GT. BEAM_CUTOFF(IBEAM) ) THEN
        WUPPER(1:2,N) = zero
        WLOWER(1:2,N) = zero
        RETURN
      ENDIF

!  set local values

      SECBAR   = AVERAGE_SECANT(N,IBEAM)
      INV_X0SQ = SECBAR * SECBAR
      MU0_INV = ONE / MU0(IBEAM)

!  Set up sum and differences for Beam source terms
!  ( sum may be required again in linearization )

      GAMMA3 = ( 2.0_dp - 3.0_dp * ASYMM_TOTAL(N) * MU0(IBEAM) ) / 4.0_dp
      GAMMA4 = One - GAMMA3

!  The reduced problem: QMAT. W = QVEC (Overwrite QVEC)
!  Restore up and down solutions

      QMAT   = EIGENVALUE(N) * EIGENVALUE(N) - INV_X0SQ
      FACTOR = OMEGA_TOTAL(n) * F1
      WVEC(1,N) = FACTOR * ( (GAMMA1(N) + MU0_INV ) * GAMMA4 + GAMMA3 * GAMMA2(N) ) / QMAT
      WVEC(2,N) = FACTOR * ( (GAMMA1(N) - MU0_INV ) * GAMMA3 + GAMMA4 * GAMMA2(N) ) / QMAT

!  Values at the layer boundaries
!  (transmittance factors have been determined in SETUPS module)

      TRANS1 = INITIAL_TRANS(N,IBEAM)
      TRANS2 = T_DELT_MUBAR(N,IBEAM) * TRANS1
      WUPPER(1:2,N) = WVEC(1:2,N)*TRANS1
      WLOWER(1:2,N) = WVEC(1:2,N)*TRANS2

!  Validates against  Toon code, 6/29/24
!write(*,*)n,GAMMA1(N),GAMMA2(N),GAMMA3,GAMMA4,WVEC(1:2,N)/FACTOR
!write(*,*)n, WUPPER(2,N), WUPPER(1,N) 
!write(*,*)n, WLOWER(2,N), WLOWER(1,N)

!  Finish

      RETURN
END SUBROUTINE EDDINGTON_BEAM_SOLUTION

end module twostream_solutions_m
