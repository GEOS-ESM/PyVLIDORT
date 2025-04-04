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
! # Regular BVP: Subroutines in this Module                     #
! #                                                             #
! #            TWOSTREAM_BVP_MATSETUP_PENTADIAG                 #
! #            TWOSTREAM_BVP_SOLUTION_PENTADIAG                 #
! #                                                             #
! #            TWOSTREAM_BVP_MATSETUP_TRIDIAG                   #
! #            TWOSTREAM_BVP_SOLUTION_TRIDIAG                   #
! #                                                             #
! ###############################################################
 
module twostream_bvproblem_m

PUBLIC

contains

SUBROUTINE TWOSTREAM_BVP_MATSETUP_PENTADIAG &
         ( MAXLAYERS, MAXTOTAL, DO_INVERSE, DO_FULLQUAD,            & ! Dimensions, flags
           DO_INCLUDE_SURFACE, DO_BRDF_SURFACE, FOURIER,            & ! input
           NLAYERS, NTOTAL, FF, SURFACE_FACTOR, STREAM_VALUE,       & ! input
           ALBEDO, BRDF_F, XPOS, XNEG, T_DELT_EIGEN, CONSSCAT, TGM, & ! input
           H_HOMP, H_HOMM, MAT, ELM, SELM,                          & ! output
           STATUS, MESSAGE )                                          ! output

      implicit none

!  precision and parameters

      INTEGER      , PARAMETER :: dp   = KIND( 1.0D0 )
      REAL(kind=dp), parameter :: zero = 0.0_dp, one = 1.0_dp, half = 0.5_dp

!  input
!  -----

!  Dimensions

      INTEGER, INTENT(IN)        :: MAXLAYERS, MAXTOTAL

!   Full quadrature control (true ==> Mu_bar = 1/sqrt(3) )

      LOGICAL      , INTENT(IN)  :: DO_FULLQUAD

!  Inverse control

      LOGICAL      , INTENT(IN)  :: DO_INVERSE
      REAL(kind=dp), INTENT(IN)  :: FF

!  Surface control

      LOGICAL, INTENT(IN)        :: DO_INCLUDE_SURFACE
      LOGICAL      , INTENT(IN)  :: DO_BRDF_SURFACE

!  control numbers

      INTEGER, INTENT(IN)        :: FOURIER
      INTEGER, INTENT(IN)        :: NLAYERS, NTOTAL

!  Surface inputs

      REAL(kind=dp), INTENT(IN)  :: SURFACE_FACTOR
      REAL(kind=dp), INTENT(IN)  :: STREAM_VALUE
      REAL(kind=dp), INTENT(IN)  :: ALBEDO
      REAL(kind=dp), INTENT(IN)  :: BRDF_F(0:1)

!  Eigenvector solutions

      REAL(kind=dp), INTENT(IN)  :: XPOS(2,MAXLAYERS)
      REAL(kind=dp), INTENT(IN)  :: XNEG(2,MAXLAYERS)

!  transmittance factors for +/- eigenvalues

      REAL(kind=dp), INTENT(IN)  :: T_DELT_EIGEN(MAXLAYERS)

!  Conservative scattering flags and parameters

      LOGICAL      , INTENT(IN)  :: CONSSCAT(MAXLAYERS)
      REAL(kind=dp), INTENT(IN)  :: TGM(MAXLAYERS)

!  Output
!  ------

!  Downwelling BOA solutions

      REAL(kind=dp), INTENT(OUT) :: H_HOMP
      REAL(kind=dp), INTENT(OUT) :: H_HOMM

!  Pentadiagonal Matrix entries for solving BCs

      REAL(kind=dp), INTENT(OUT) :: MAT(MAXTOTAL,5)

!  Pentadiagonal elimination marix

      REAL(kind=dp), INTENT(OUT) :: ELM (MAXTOTAL,4)

!  single layer elimination matrix

      REAL(kind=dp), INTENT(OUT) :: SELM (2,2)

!  status

      INTEGER, INTENT(OUT)        :: STATUS
      CHARACTER*(*), INTENT(OUT)  :: MESSAGE

!  local variables
!  ---------------

!  square matrix for the single layer case

      REAL(kind=dp)     :: SMAT (2,2)

!  Help

      INTEGER           :: I, N, N1, NM, NP, INM, INP
      REAL(kind=dp)     :: XPNET, XMNET, FACTOR, BET, DEN, R2_HOMP, R2_HOMM, MFAC
      CHARACTER(LEN=3)  :: CI

!  Stability check value

      REAL(kind=dp)     :: SMALLNUM=1.0D-20

!  Status

      status = 0
      message = ' '

!  Additional setups for the lowest layer
!  For Lambertian reflectance, all streams are the same
!  For BRDF, code added 4 May 2009 by R. Spurr

!  Rob Fix 2p4a, 5/26/17. Be careful with energy conservation

      MFAC = SURFACE_FACTOR * HALF
!      MFAC = SURFACE_FACTOR * STREAM_VALUE
!      if ( DO_FULLQUAD ) MFAC = HALF * SURFACE_FACTOR

      if ( CONSSCAT(NLAYERS) ) then
         H_HOMP = ONE - TGM(NLAYERS)
         H_HOMM = TGM(NLAYERS)
      else
         H_HOMP = XPOS(1,NLAYERS)
         H_HOMM = XNEG(1,NLAYERS)
      endif

      R2_HOMP = zero ; R2_HOMM = zero
      IF ( DO_INCLUDE_SURFACE ) THEN
        IF ( DO_BRDF_SURFACE ) THEN
          FACTOR = MFAC * BRDF_F(FOURIER)
        ELSE
          FACTOR = MFAC * ALBEDO
        ENDIF
        R2_HOMP = FACTOR * H_HOMP
        R2_HOMM = FACTOR * H_HOMM
      ENDIF

!  Inclusion of surface contribution in BV Problem matrix

      if ( CONSSCAT(NLAYERS) ) then
         XPNET = - TGM(NLAYERS)
         XMNET = ONE + TGM(NLAYERS)
      else
         XPNET = XPOS(2,NLAYERS)
         XMNET = XNEG(2,NLAYERS)
      endif

      IF ( DO_INCLUDE_SURFACE ) THEN
        XPNET = XPNET - R2_HOMP 
        XMNET = XMNET - R2_HOMM
      ENDIF

!  set up BVP matrix
!  -----------------

!  If Nlayers = 1, special case, ordinary 2x2

      IF ( NLAYERS .EQ. 1 ) THEN

!  Only Top and bottom BC (with surface reflection)

        if ( CONSSCAT(NLAYERS) ) then
          SMAT(1,1) = one            !  downwelling at top of slab (x=0).   alpha
          SMAT(1,2) = zero           !  downwelling at top of slab (x=0).   beta
          SMAT(2,1) = XPNET          !  upwelling   at bot of slab (x=tau). alpha
          SMAT(2,2) = XMNET          !  upwelling   at bot of slab (x=tau). bet
        else
          SMAT(1,1) = XPOS(1,NLAYERS)
          SMAT(1,2) = XNEG(1,NLAYERS) * T_DELT_EIGEN(NLAYERS)
          SMAT(2,1) = XPNET * T_DELT_EIGEN(NLAYERS)
          SMAT(2,2) = XMNET
        endif

! write(*,*)SMAT(1,1:2),SMAT(2,1:2)

!  If NLAYERS > 1, set up Pentadiagonal matrix

      ELSE

!  Zero for both Fourier components (Important bug!)

        MAT = zero
        ELM = zero

!  top BC for layer 1: no downward diffuse radiation
!  intermediate layer boundaries
!  bottom BC (including surface reflected term)

        if ( .NOT. do_inverse ) then
          if ( CONSSCAT(1) ) then
            MAT(1,3)  = one                               !  downwelling at top of layer (x=0).   alpha
            MAT(1,4)  = zero                              !  downwelling at top of layer (x=0).   beta
          else
            MAT(1,3)  = XPOS(1,1)                         !  downwelling at top of layer (x=0).   alpha
            MAT(1,4)  = XNEG(1,1) * T_DELT_EIGEN(1)       !  downwelling at top of layer (x=0).   beta
          endif
          DO N = 2, NLAYERS
            N1 =  N - 1 ; NM = 2*N1 ; NP = NM + 1
            IF ( CONSSCAT(N1) ) THEN
              MAT(NM,2) =   one - tgm(n1)                 !   dnwelling at bottom of layer N1.   alpha
              MAT(NM,3) =   tgm(n1)                       !   dnwelling at bottom of layer N1.   beta
              MAT(NP,1) =   - tgm(n1)                     !   Upwelling at bottom of layer N1.   alpha
              MAT(NP,2) =   one + tgm(n1)                 !   Upwelling at bottom of layer N1.   beta
            ELSE
              MAT(NM,2) =   XPOS(1,N1) * T_DELT_EIGEN(N1) !   dnwelling at bottom of layer N1.   alpha
              MAT(NM,3) =   XNEG(1,N1)                    !   dnwelling at bottom of layer N1.   beta
              MAT(NP,1) =   XPOS(2,N1) * T_DELT_EIGEN(N1) !   Upwelling at bottom of layer N1.   alpha
              MAT(NP,2) =   XNEG(2,N1)                    !   Upwelling at bottom of layer N1.   beta
            ENDIF
            IF ( CONSSCAT(N) ) THEN
              MAT(NM,4) =  - one                          ! - dnwelling at Top    of layer N.    alpha
              MAT(NM,5) = zero                            ! - dnwelling at Top    of layer N.    beta
              MAT(NP,3) = zero                            ! - Upwelling at Top    of layer N.    alpha
              MAT(NP,4) = - one                           ! - Upwelling at Top    of layer N.    beta
            ELSE
              MAT(NM,4) = - XPOS(1,N)                     ! - dnwelling at Top    of layer N.    alpha
              MAT(NM,5) = - XNEG(1,N)  * T_DELT_EIGEN(N)  ! - dnwelling at Top    of layer N.    beta
              MAT(NP,3) = - XPOS(2,N)                     ! - Upwelling at Top    of layer N.    alpha
              MAT(NP,4) = - XNEG(2,N)  * T_DELT_EIGEN(N)  ! - Upwelling at Top    of layer N.    beta
            ENDIF
          ENDDO
          if ( CONSSCAT(NLAYERS) ) then
            MAT(NTOTAL,2) = XPNET  !  upwelling at surface. alpha
            MAT(NTOTAL,3) = XMNET  !  upwelling at surface. beta
          else
            MAT(NTOTAL,2) = XPNET * T_DELT_EIGEN(NLAYERS) !  upwelling at surface. alpha
            MAT(NTOTAL,3) = XMNET                         !  upwelling at surface. beta
          endif
        endif

!  Inverted set

        if ( do_inverse ) then
          if ( CONSSCAT(NLAYERS) ) then
            MAT(1,3) = XMNET                                !  upwelling at surface. beta
            MAT(1,4) = XPNET                                !  upwelling at surface. alpha
          else
            MAT(1,3)  = XMNET                               !  upwelling at surface. beta
            MAT(1,4)  = XPNET * T_DELT_EIGEN(NLAYERS)       !  upwelling at surface. alpha
          endif
          DO N = 2, NLAYERS
            N1 =  N - 1 ;  INM = NTOTAL - 2*N1  ; INP = INM + 1
            IF ( CONSSCAT(N1) ) THEN
              MAT(INM,4) =   tgm(n1)                        !   Upwelling at bottom of layer N1.   beta
              MAT(INM,5) =   one - tgm(n1)                  !   Upwelling at bottom of layer N1.   alpha
              MAT(INP,3) =   one + tgm(n1)                  !   dnwelling at bottom of layer N1.   beta
              MAT(INP,4) =   - tgm(n1)                      !   dnwelling at bottom of layer N1.   alpha
            ELSE
              MAT(INM,4) =    XNEG(2,N1)                    !   Upwelling at bottom of layer N1.   beta
              MAT(INM,5) =    XPOS(2,N1) * T_DELT_EIGEN(N1) !   Upwelling at bottom of layer N1.   alpha
              MAT(INP,3) =    XNEG(1,N1)                    !   dnwelling at bottom of layer N1.   beta
              MAT(INP,4) =    XPOS(1,N1) * T_DELT_EIGEN(N1) !   dnwelling at bottom of layer N1.   alpha
            ENDIF
            IF ( CONSSCAT(N) ) THEN
              MAT(INM,2) =   zero                           ! - Upwelling at Top    of layer N.    beta
              MAT(INM,3) = - one                            ! - Upwelling at Top    of layer N.    alpha
              MAT(INP,1) = - one                            ! - dnwelling at Top    of layer N.    beta
              MAT(INP,2) =   zero                           ! - dnwelling at Top    of layer N.    alpha
            ELSE
              MAT(INM,2) =  - XNEG(2,N)  * T_DELT_EIGEN(N)  ! - Upwelling at Top    of layer N.    beta
              MAT(INM,3) =  - XPOS(2,N)                     ! - Upwelling at Top    of layer N.    alpha
              MAT(INP,1) =  - XNEG(1,N)  * T_DELT_EIGEN(N)  ! - dnwelling at Top    of layer N.    beta
              MAT(INP,2) =  - XPOS(1,N)                     ! - dnwelling at Top    of layer N.    alpha
            ENDIF
          ENDDO
          if ( CONSSCAT(1) ) then
            MAT(NTOTAL,2)  = zero                           !  downwelling at top of layer (x=0).   beta
            MAT(NTOTAL,3)  = one                            !  downwelling at top of layer (x=0).   alpha
          else
            MAT(NTOTAL,2) = XNEG(1,1) * T_DELT_EIGEN(1)     !  downwelling at top of layer (x=0).   beta
            MAT(NTOTAL,3) = XPOS(1,1)                       !  downwelling at top of layer (x=0).   alpha
          endif
        endif

      ENDIF

!  Scaling

      MAT = FF * MAT

!  Elimination of BVP pentadiagonal matrix
!  ---------------------------------------

      IF ( NLAYERS .GT. 1 ) THEN

!  Row 1

        ELM(1,4) = zero
        ELM(1,3) = one / MAT(1,3)
        ELM(1,1) = - MAT(1,4) * ELM(1,3)
        ELM(1,2) = - MAT(1,5) * ELM(1,3)

!  Row 2; includes first check for singularity

        ELM(2,4) = zero
        bet = MAT(2,3) + MAT(2,2)*ELM(1,1)
!write(*,*)FOURIER,MAT(2,3), MAT(2,2)*ELM(1,1),bet
        IF ( DABS(Bet) .LT. SMALLNUM ) THEN
          message = 'Singularity in Pentadiagonal Matrix, Row #  2'
          status = 1
          return
        endif
        bet = - one / bet
        ELM(2,1) = (MAT(2,4) + MAT(2,2)*ELM(1,2)) * bet
        ELM(2,2) = MAT(2,5) * bet
        ELM(2,3) = bet
        ELM(2,4) = zero

!  Rows 3-NT: reduce to upper triangular; includes checks for singularity

        do i = 3, ntotal
          bet = MAT(i,2) + MAT(i,1) * ELM(i-2,1)
          den = MAT(i,3) + MAT(i,1) * ELM(i-2,2) + bet * ELM(i-1,1)
          IF ( DABS(DEN) .LT. SMALLNUM ) THEN
            WRITE(CI, '(I3)' ) I
            message = 'Singularity in Pentadiagonal Matrix, Row #'//CI
            status = 1
            return
          endif
          den = - one / den
          ELM(i,1) = (MAT(i,4) + bet*ELM(i-1,2)) * den
          ELM(i,2) = MAT(i,5) * den
          ELM(i,3) = bet
          ELM(i,4) = den
        enddo

!  Elimination for Single Layer only
!  ----------------------------------

      ELSE IF ( NLAYERS .EQ. 1 ) THEN

        DEN = SMAT(1,1)*SMAT(2,2) - SMAT(1,2)*SMAT(2,1)
        IF ( DABS(DEN) .LT. SMALLNUM ) THEN
          message = 'Singularity in 1-layer 2x2 Matrix'
          status = 1
          return
        ENDIF
        DEN = one / DEN
        SELM(1,1) =   SMAT(2,2) * DEN
        SELM(1,2) = - SMAT(1,2) * DEN
        SELM(2,1) = - SMAT(2,1) * DEN
        SELM(2,2) =   SMAT(1,1) * DEN

      ENDIF

!  finish

      RETURN
 END SUBROUTINE TWOSTREAM_BVP_MATSETUP_PENTADIAG

!

SUBROUTINE TWOSTREAM_BVP_SOLUTION_PENTADIAG &
      ( MAXLAYERS, MAXBEAMS, MAXTOTAL,                   & ! Dimensions
        DO_INVERSE, DO_FULLQUAD,                         & ! Flags
        DO_INCLUDE_SURFACE, DO_INCLUDE_DIRECTBEAM,       & ! inputs
        DO_INCLUDE_SURFEMISS, DO_BRDF_SURFACE,           & ! inputs
        FOURIER, IPARTIC, NLAYERS, NTOTAL,               & ! inputs
        FF, SURFACE_FACTOR, STREAM_VALUE,                & ! inputs
        ALBEDO, BRDF_F, EMISS, SURFBB,                   & ! inputs
        DIRECT_BEAM, XPOS, XNEG, WUPPER, WLOWER,         & ! inputs
        MAT, ELM, SELM, CONSSCAT,                        & ! inputs
        H_PARTIC, LCON, MCON, LCON_XVEC, MCON_XVEC  )      ! Output

      implicit none

!  precision and parameters

      INTEGER      , PARAMETER :: dp   = KIND( 1.0D0 )
      REAL(kind=dp), parameter :: zero = 0.0_dp, half = 0.5_dp

!  Input arguments
!  ---------------

!  Dimensions

      INTEGER, INTENT(IN)        :: MAXBEAMS, MAXLAYERS, MAXTOTAL

!   Full quadrature control (true ==> Mu_bar = 1/sqrt(3) )

      LOGICAL      , INTENT(IN)  :: DO_FULLQUAD

!  Inverse control

      LOGICAL      , INTENT(IN)  :: DO_INVERSE
      REAL(kind=dp), INTENT(IN)  :: FF

!  Inclusion flags for surface terms

      LOGICAL, INTENT(IN)        :: DO_INCLUDE_DIRECTBEAM
      LOGICAL, INTENT(IN)        :: DO_INCLUDE_SURFACE
      LOGICAL, INTENT(IN)        :: DO_BRDF_SURFACE
      LOGICAL, INTENT(IN)        :: DO_INCLUDE_SURFEMISS

!  Surface inputs

      REAL(kind=dp), INTENT(IN)  :: SURFACE_FACTOR
      REAL(kind=dp), INTENT(IN)  :: STREAM_VALUE
      REAL(kind=dp), INTENT(IN)  :: ALBEDO
      REAL(kind=dp), INTENT(IN)  :: BRDF_F(0:1)
      REAL(kind=dp), INTENT(IN)  :: EMISS, SURFBB

!  Fourier component and beam number, layers

      INTEGER, INTENT(IN)        :: FOURIER, IPARTIC
      INTEGER, INTENT(IN)        :: NLAYERS, NTOTAL

!  Direct beam

      REAL(kind=dp), INTENT(IN)  :: DIRECT_BEAM ( MAXBEAMS )

!  Eigenvector solutions

      REAL(kind=dp), INTENT(IN)  :: XPOS(2,MAXLAYERS)
      REAL(kind=dp), INTENT(IN)  :: XNEG(2,MAXLAYERS)

!  Particular solutions

      REAL(kind=dp), INTENT(IN)  :: WLOWER ( 2, MAXLAYERS )
      REAL(kind=dp), INTENT(IN)  :: WUPPER ( 2, MAXLAYERS )

!  Pentadiagonal Matrix entries for solving BCs

      REAL(kind=dp), INTENT(IN)  :: MAT(MAXTOTAL,5)

!  Pentadiagonal elimination matrix

      REAL(kind=dp), INTENT(IN)  :: ELM (MAXTOTAL,4)

!  Single layer elimination matrix

      REAL(kind=dp), INTENT(IN)  :: SELM (2,2)

!  Conservative scattering flags

      LOGICAL      , INTENT(IN)  :: CONSSCAT(MAXLAYERS)

!  Output
!  ------

!  Downwelling BOA solution

      REAL(kind=dp), INTENT(OUT) :: H_PARTIC

!  Solution constants of integration, and related quantities

      REAL(kind=dp), INTENT(OUT) :: LCON(MAXLAYERS)
      REAL(kind=dp), INTENT(OUT) :: MCON(MAXLAYERS)

      REAL(kind=dp), INTENT(OUT) :: LCON_XVEC(2,MAXLAYERS)
      REAL(kind=dp), INTENT(OUT) :: MCON_XVEC(2,MAXLAYERS)

!  Local variables
!  ---------------

!  Column vectors for solving BCs. Not saved.

      REAL(kind=dp)       :: COL    (MAXTOTAL)
      REAL(kind=dp)       :: SCOL   (2)

!  Other variables

      INTEGER             :: N, N1, I, NM, NP, INP, INM, NI, NLAY1
      REAL(kind=dp)       :: FACTOR, DEN, R2_PARTIC, TOA_TERM, SURFACE_TERM, MFAC
      REAL(kind=dp)       :: NEW_SCOL1

!  Additional setups for the surface and TOA levels
!  ------------------------------------------------

!  Zero total reflected contribution (R2_PARTIC) before calculation
!  For Lambertian reflectance, all streams are the same
!  For BRDF, code added 4 May 2009 by R. Spurr

!  Rob Fix 2p4a, 5/26/17. Be careful with energy conservation

      MFAC = SURFACE_FACTOR * HALF
!      MFAC = SURFACE_FACTOR * STREAM_VALUE
!      if ( DO_FULLQUAD ) MFAC = HALF * SURFACE_FACTOR

      R2_PARTIC = zero
      H_PARTIC = WLOWER(1,NLAYERS)

      IF ( DO_INCLUDE_SURFACE ) THEN
        IF ( DO_BRDF_SURFACE ) THEN
          FACTOR = MFAC * BRDF_F(FOURIER)
        ELSE
          FACTOR = MFAC * ALBEDO
        ENDIF
        R2_PARTIC = H_PARTIC * FACTOR
      ENDIF

!  Surface Term------------
!    lowest (surface) boundary with albedo (diffuse radiation terms only)
!    with non-zero albedo, include integrated downward reflectances
!    no albedo, similar code excluding integrated reflectance
!    Add direct beam solution (only to final level)
!    Add thermal emission of ground surface (only to final level)
!mick fix 2/14/2012 - changed treatment of emissivity to be consistent
!                     with LIDORT & VLIDORT

      SURFACE_TERM = - WLOWER(2,NLAYERS)
      IF ( DO_INCLUDE_SURFACE ) THEN
        SURFACE_TERM = SURFACE_TERM + R2_PARTIC
        IF ( DO_INCLUDE_DIRECTBEAM ) THEN
          SURFACE_TERM = SURFACE_TERM + DIRECT_BEAM(IPARTIC)
        ENDIF
      ENDIF
      IF ( DO_INCLUDE_SURFEMISS ) THEN
        SURFACE_TERM = SURFACE_TERM + SURFBB * EMISS
      ENDIF
!  Upper boundary for layer 1: no downward diffuse radiation

      TOA_TERM  = - WUPPER(1,1)

!  --set up Column for solution vector (the "B" as in AX=B)
!  --------------------------------------------------------

!  Pentadiagonal case

      IF ( NLAYERS .GT. 1 ) THEN

!  zero column vector

        COL(1:NTOTAL) = zero

!  Fill vector

        if ( do_inverse ) then
          COL(1)  = SURFACE_TERM
          DO N = 2, NLAYERS
            N1 = N - 1
            INM = NTOTAL - 2*N1 ; INP = INM + 1
            COL(INP) = WUPPER(1,N) - WLOWER(1,N1)
            COL(INM) = WUPPER(2,N) - WLOWER(2,N1)
          ENDDO
          COL(NTOTAL) = TOA_TERM
        else
          COL(1)  = TOA_TERM
          DO N = 2, NLAYERS
            N1 = N - 1
            NM = 2*N1 ; NP = NM + 1
            COL(NM) = WUPPER(1,N) - WLOWER(1,N1)
            COL(NP) = WUPPER(2,N) - WLOWER(2,N1)
          ENDDO
          COL(NTOTAL) = SURFACE_TERM
        endif

!  Scaling

!mick fix 5/29/2015 - tailor dimensions
        !COL = FF * COL
        COL(1:NTOTAL) = FF * COL(1:NTOTAL)

!  debug
!        do n = 1, ntotal
!          write(44,'(2i3,1p3e24.12)')n,1,COL(N)
!        enddo
!       pause

!  If Nlayers = 1, special case

      ELSE IF ( NLAYERS .EQ. 1 ) THEN
        SCOL(1) = TOA_TERM
        SCOL(2) = SURFACE_TERM
      ENDIF

!  --Solve the boundary problem for this Fourier component (back substitution)
! ---------------------------------------------------------------------------

!  Pentadiagonal back-substitution

      IF ( NLAYERS .GT. 1 ) THEN

!  Fill up back-substitution array

        COL(1) = COL(1) * ELM(1,3)
        COL(2) = (MAT(2,2)*COL(1) - COL(2)) * ELM(2,3)
        do I = 3, NTOTAL
          DEN = ELM(i,4)
          COL(i) = (MAT(i,1)*COL(i-2)+ELM(i,3)*COL(i-1)-COL(i)) * DEN
        enddo

!  Back-substitution

        N1 = NTOTAL-1
        COL(N1) = COL(N1) + ELM(N1,1) * COL(NTOTAL)
        do i = NTOTAL-2, 1, -1
          COL(i) = COL(i) + ELM(i,1) * COL(i+1) + ELM(i,2) * COL(i+2)
        enddo

!  Set integration constants LCON and MCON for -/+ eigensolutions, all layers

        if ( do_inverse ) then
           NLAY1 = 1 + NLAYERS
           DO N = 1, NLAYERS
              NI = NLAY1 - N
              INP = 2*NI ; INM = INP - 1
              LCON(N) = COL(INP)
              MCON(N) = COL(INM)
           ENDDO
        else
           DO N = 1, NLAYERS
              NM = 2*N-1 ; NP = NM + 1
              LCON(N) = COL(NM)
              MCON(N) = COL(NP)
           ENDDO
        endif

!  Solve the boundary problem: Single Layer only
!mick fix 1/9/2011 - defined NEW_SCOL1 so a MODIFIED version of SCOL(1)
!                    is not being used in computing SCOL(2)

      ELSE IF ( NLAYERS .EQ. 1 ) THEN
        NEW_SCOL1 = SELM(1,1) * SCOL(1) + SELM(1,2) * SCOL(2)
        SCOL(2)   = SELM(2,1) * SCOL(1) + SELM(2,2) * SCOL(2)
        LCON(1) = NEW_SCOL1
        MCON(1) = SCOL(2)
      ENDIF

!  debug
!      if ( fourier.eq.0 .and.ipartic.eq.1) then
!        do n = 1, nlayers
!          write(34,'(i3,1p2e24.12)')n,LCON(N),MCON(N)
!        enddo
!      endif

!  Associated quantities
!  ---------------------

      DO N = 1, NLAYERS
        IF ( .not. CONSSCAT(n) ) then
          LCON_XVEC(1:2,N) = LCON(N)*XPOS(1:2,N)
          MCON_XVEC(1:2,N) = MCON(N)*XNEG(1:2,N)
        endif
      ENDDO

!  Finish

      RETURN
END SUBROUTINE TWOSTREAM_BVP_SOLUTION_PENTADIAG

!

SUBROUTINE TWOSTREAM_BVP_MATSETUP_TRIDIAG &
         ( MAXLAYERS, MAXTOTAL, DO_FULLQUAD,          & ! Input
           DO_INCLUDE_SURFACE, DO_BRDF_SURFACE,       & ! Input
           FOURIER, NLAYERS, NTOTAL, SURFACE_FACTOR,  & ! Input
           STREAM_VALUE, ALBEDO, BRDF_F, XPOS,        & ! Input
           GAMMA, T_DELT_EIGEN, CONSSCAT, TGM,        & ! input
           AMAT_TD, BMAT_TD, DMAT_TD, EMAT_TD, SELM,  & ! output
           STATUS, MESSAGE )                            ! output

!  Conservative case does not yield to this TRIDIAG treatment

      implicit none

!  precision and parameters

      INTEGER      , PARAMETER :: dp   = KIND( 1.0D0 )
      REAL(kind=dp), parameter :: zero = 0.0_dp, one = 1.0_dp, half = 0.5_dp

!  input
!  -----

!  Dimensions

      INTEGER, INTENT(IN)        :: MAXLAYERS, MAXTOTAL

!   Full quadrature control (true ==> Mu_bar = 1/sqrt(3) )

      LOGICAL      , INTENT(IN)  :: DO_FULLQUAD

!  Surface control

      LOGICAL      , INTENT(IN)  :: DO_INCLUDE_SURFACE
      LOGICAL      , INTENT(IN)  :: DO_BRDF_SURFACE

!  control integers

      INTEGER, INTENT(IN)        :: FOURIER
      INTEGER, INTENT(IN)        :: NLAYERS, NTOTAL

!  Surface inputs

      REAL(kind=dp), INTENT(IN)  :: SURFACE_FACTOR
      REAL(kind=dp), INTENT(IN)  :: ALBEDO, STREAM_VALUE
      REAL(kind=dp), INTENT(IN)  :: BRDF_F(0:1)

!  Eigenvector solutions, Gamma functions
!  transmittance factors for +/- eigenvalues
!         GAMMA = XPOS(2,N) / XPOS(1,N)

      REAL(kind=dp), INTENT(IN)  :: XPOS(2,MAXLAYERS)
      REAL(kind=dp), INTENT(IN)  :: GAMMA(MAXLAYERS)
      REAL(kind=dp), INTENT(IN)  :: T_DELT_EIGEN(MAXLAYERS)

!  Conservative scattering flags and parameters.

      LOGICAL      , INTENT(IN)  :: CONSSCAT(MAXLAYERS)
      REAL(kind=dp), INTENT(IN)  :: TGM(MAXLAYERS)

!  Output
!  ------

!  Tridiagonal Matrix entries for solving BCs

      REAL(kind=dp), INTENT(OUT) :: AMAT_TD(MAXTOTAL)
      REAL(kind=dp), INTENT(OUT) :: BMAT_TD(MAXTOTAL)
      REAL(kind=dp), INTENT(OUT) :: DMAT_TD(MAXTOTAL)

!  Tridiagonal help matrices

      REAL(kind=dp), INTENT(OUT) :: EMAT_TD (MAXLAYERS,4)

!  single layer elimination matrix

      REAL(kind=dp), INTENT(OUT) :: SELM (2,2)

!  status

      INTEGER, INTENT(OUT)       :: STATUS
      CHARACTER*(*), INTENT(OUT) :: MESSAGE

!  local variables
!  ---------------

!  square matrix for the single layer case

      REAL(kind=dp)     :: SMAT (2,2)

!  Help

      INTEGER           :: J, N
      REAL(kind=dp)     :: REFLEC, MFAC, DEN !, gammad(maxlayers)

!  Stability check value

      REAL(kind=dp)     :: SMALLNUM=1.0D-20

!  Status   

      status = 0
      message = ' '

!  Initialize

      AMAT_TD = zero ; BMAT_TD = zero
      DMAT_TD = zero ; EMAT_TD = zero

!  Surface

!  Rob Fix 2p4a, 5/26/17. Be careful with energy conservation

      MFAC = HALF * SURFACE_FACTOR

      REFLEC = zero
      IF ( DO_INCLUDE_SURFACE ) THEN
        IF ( DO_BRDF_SURFACE ) THEN
          REFLEC = MFAC * BRDF_F(FOURIER)
        ELSE
          REFLEC = MFAC * ALBEDO
        ENDIF
      ENDIF

!  Set up Ematrices
!  ----------------

      DO N = 1, NLAYERS
        IF ( CONSSCAT(N) ) THEN
          EMAT_TD(n,1) = ONE + TGM(N)
          EMAT_TD(n,2) = - TGM(N)
          EMAT_TD(n,3) = TGM(N)
          EMAT_TD(n,4) = ONE - TGM(N)
        ELSE
          EMAT_TD(n,1) = ONE + GAMMA(N) * T_DELT_EIGEN(N)
          EMAT_TD(n,2) = ONE - GAMMA(N) * T_DELT_EIGEN(N)
          EMAT_TD(n,3) = GAMMA(N) + T_DELT_EIGEN(N)
          EMAT_TD(n,4) = GAMMA(N) - T_DELT_EIGEN(N)
        ENDIF
      ENDDO

!  Check out #14
!DO N = 1, NLAYERS
!  write(*,*)n,emat_td(n,1:4)
!ENDDO

! open(44,file='GammaCheat.dat',status='old')
! do n = 1, nlayers
!    read(44,*)j,mfac,gammad(n),mfac
!          EMAT_TD(n,1) = ONE + GAMMAd(N) * T_DELT_EIGEN(N)
!          EMAT_TD(n,2) = ONE - GAMMAd(N) * T_DELT_EIGEN(N)
!          EMAT_TD(n,3) = GAMMAd(N) + T_DELT_EIGEN(N)
!          EMAT_TD(n,4) = GAMMAd(N) - T_DELT_EIGEN(N)
!enddo
!close(44)

!  If Nlayers = 1, special case, ordinary 2x2

      IF ( NLAYERS .EQ. 1 ) THEN

!  Only Top and bottom BC (with surface reflection)

        IF ( CONSSCAT(NLAYERS) ) then

          SMAT(1,1) = zero
          SMAT(1,2) = one
          SMAT(2,1) = emat_td(nlayers,1) - REFLEC * emat_td(nlayers,3)
          SMAT(2,2) = emat_td(nlayers,2) - REFLEC * emat_td(nlayers,4)

        ELSE

          SMAT(1,1) = emat_td(nlayers,1)
          SMAT(1,2) = - emat_td(nlayers,2)
          SMAT(2,1) = emat_td(nlayers,1) - REFLEC * emat_td(nlayers,3)
          SMAT(2,2) = emat_td(nlayers,2) - REFLEC * emat_td(nlayers,4)

        ENDIF

        DEN = SMAT(1,1)*SMAT(2,2) - SMAT(1,2)*SMAT(2,1)
        IF ( DABS(DEN) .LT. SMALLNUM ) THEN
          message = 'Singularity in 1-layer 2x2 Matrix'
          status = 1
          return
        ENDIF
        DEN = one / DEN
        SELM(1,1) =   SMAT(2,2) * DEN
        SELM(1,2) = - SMAT(1,2) * DEN
        SELM(2,1) = - SMAT(2,1) * DEN
        SELM(2,2) =   SMAT(1,1) * DEN

!  If NLAYERS > 1, set up tridiagonal matrix

      ELSE

!  Set up A,B,D matrices
!  ---------------------

      IF ( CONSSCAT(1) ) THEN
        AMAT_TD(1) = zero
        BMAT_TD(1) = zero
        DMAT_TD(1) = one
      ELSE
        AMAT_TD(1) = zero
        BMAT_TD(1) = EMAT_TD(1,1)
        DMAT_TD(1) = - EMAT_TD(1,2)
      ENDIF

      DO N = 1, NLAYERS-1
      IF ( CONSSCAT(N+1) ) THEN
           j = 2*n+1
           AMAT_TD(j) = emat_td(n,2)*emat_td(n,3)   - emat_td(n,4)*emat_td(n,1)
           BMAT_TD(j) = - emat_td(n,3)
           DMAT_TD(j) = emat_td(n,1)
           j = 2*n
           AMAT_TD(j) = emat_td(n,1)
           BMAT_TD(j) = emat_td(n,2)
           DMAT_TD(j) = - one
      ELSE
           j = 2*n+1
           AMAT_TD(j) = emat_td(n,2)*emat_td(n,3)   - emat_td(n,4)*emat_td(n,1)
           BMAT_TD(j) = emat_td(n,1)*emat_td(n+1,1) - emat_td(n,3)*emat_td(n+1,3)
           DMAT_TD(j) = emat_td(n,3)*emat_td(n+1,4) - emat_td(n,1)*emat_td(n+1,2)
           j = 2*n
           AMAT_TD(j) = emat_td(n+1,2)*emat_td(n,1)   - emat_td(n,3)  *emat_td(n+1,4)
           BMAT_TD(j) = emat_td(n,2)  *emat_td(n+1,2) - emat_td(n,4)  *emat_td(n+1,4)
           DMAT_TD(j) = emat_td(n+1,1)*emat_td(n+1,4) - emat_td(n+1,2)*emat_td(n+1,3)
      ENDIF
      ENDDO
      j = NTOTAL
      AMAT_TD(j) = emat_td(nlayers,1) - REFLEC * emat_td(nlayers,3)
      BMAT_TD(j) = emat_td(nlayers,2) - REFLEC * emat_td(nlayers,4)
      DMAT_TD(j) = zero

!write(*,*)'check # 16'
!do j = 1, ntotal
!   write(*,*)j,amat_td(j),BMAT_TD(j),DMAT_TD(j)
!enddo

        ENDIF
!  finish

      RETURN
END SUBROUTINE TWOSTREAM_BVP_MATSETUP_TRIDIAG

!

SUBROUTINE TWOSTREAM_BVP_SOLUTION_TRIDIAG &
      ( MAXLAYERS, MAXBEAMS, MAXTOTAL, DO_FULLQUAD,        & ! Dimensions
        DO_INCLUDE_SURFACE, DO_INCLUDE_DIRECTBEAM,         & ! inputs
        DO_INCLUDE_SURFEMISS, DO_BRDF_SURFACE,             & ! inputs
        FOURIER, IPARTIC, NLAYERS, NTOTAL, SURFACE_FACTOR, & ! inputs
        STREAM_VALUE, ALBEDO, BRDF_F, EMISS, SURFBB,       & ! inputs
        DIRECT_BEAM, XPOS, XNEG, WUPPER, WLOWER, CONSSCAT, & ! inputs
        GAMMA, T_DELT_EIGEN, AMAT_TD, BMAT_TD, DMAT_TD, EMAT_TD, SELM, & ! inputs
        LCON, MCON, LCON_XVEC, MCON_XVEC  )                  ! Output

      implicit none

!  precision and parameters

      INTEGER      , PARAMETER :: dp   = KIND( 1.0D0 )
      REAL(kind=dp), parameter :: zero = 0.0_dp, one  = 1.0_dp, half = 0.5_dp

!  Input arguments
!  ---------------

!  Dimensions

      INTEGER, INTENT(IN)        :: MAXBEAMS, MAXLAYERS, MAXTOTAL

!   Full quadrature control (true ==> Mu_bar = 1/sqrt(3) )

      LOGICAL      , INTENT(IN)  :: DO_FULLQUAD

!  Inclusion flags for Surface terms

      LOGICAL, INTENT(IN)        :: DO_INCLUDE_DIRECTBEAM
      LOGICAL, INTENT(IN)        :: DO_INCLUDE_SURFACE
      LOGICAL, INTENT(IN)        :: DO_BRDF_SURFACE
      LOGICAL, INTENT(IN)        :: DO_INCLUDE_SURFEMISS

!  Surface inputs

      REAL(kind=dp), INTENT(IN)  :: SURFACE_FACTOR
      REAL(kind=dp), INTENT(IN)  :: ALBEDO, STREAM_VALUE
      REAL(kind=dp), INTENT(IN)  :: BRDF_F(0:1)
      REAL(kind=dp), INTENT(IN)  :: EMISS, SURFBB

!  Fourier component and beam number, layers

      INTEGER, INTENT(IN)        :: FOURIER, IPARTIC
      INTEGER, INTENT(IN)        :: NLAYERS, NTOTAL

!  Direct beam

      REAL(kind=dp), INTENT(IN)  :: DIRECT_BEAM ( MAXBEAMS )

!  Eigenvector solutions

      REAL(kind=dp), INTENT(IN)  :: XPOS(2,MAXLAYERS)
      REAL(kind=dp), INTENT(IN)  :: XNEG(2,MAXLAYERS)

!  Gamma

      REAL(kind=dp), INTENT(IN)  :: GAMMA(MAXLAYERS)

!  transmittance factors for +/- eigenvalues

      REAL(kind=dp), INTENT(IN)  :: T_DELT_EIGEN(MAXLAYERS)

!  Particular solutions

      REAL(kind=dp), INTENT(IN)  :: WLOWER ( 2, MAXLAYERS )
      REAL(kind=dp), INTENT(IN)  :: WUPPER ( 2, MAXLAYERS )


!  Tridiagonal Matrix entries for solving BCs

      REAL(kind=dp), INTENT(IN) :: AMAT_TD(MAXTOTAL)
      REAL(kind=dp), INTENT(IN) :: BMAT_TD(MAXTOTAL)
      REAL(kind=dp), INTENT(IN) :: DMAT_TD(MAXTOTAL)

!  Tridiagonal help matrices

      REAL(kind=dp), INTENT(IN) :: EMAT_TD (MAXLAYERS,4)

!  Single layer elimination matrix

      REAL(kind=dp), INTENT(IN) :: SELM (2,2)

!  Conservative scattering flags. 

      LOGICAL      , INTENT(IN)  :: CONSSCAT(MAXLAYERS)

!  Output
!  ------

!  Solution constants of integration, and related quantities

      REAL(kind=dp), INTENT(OUT) :: LCON(MAXLAYERS)
      REAL(kind=dp), INTENT(OUT) :: MCON(MAXLAYERS)

      REAL(kind=dp), INTENT(OUT) :: LCON_XVEC(2,MAXLAYERS)
      REAL(kind=dp), INTENT(OUT) :: MCON_XVEC(2,MAXLAYERS)

!  Local variables
!  ---------------

!  Column vectors for solving BCs. Not saved.

      REAL(kind=dp)       :: COL_TD    (MAXTOTAL)

!  Auxiliary matrices

      REAL(kind=dp)       :: AS(MAXTOTAL)
      REAL(kind=dp)       :: DS(MAXTOTAL)
      REAL(kind=dp)       :: X(MAXTOTAL)

!  Tridiagonal solutions

      REAL(kind=dp)       :: Y(MAXTOTAL)
      REAL(kind=dp)       :: Y1(MAXLAYERS)
      REAL(kind=dp)       :: Y2(MAXLAYERS)

!  Other variables

      INTEGER             :: N, J, J1, J2, I
      REAL(kind=dp)       :: REFLEC, SURFACE_TERM, DEN, K1, K2, MFAC !, GAMMAd(maxlayers)

!  Additional setups for the surface and TOA levels
!  ------------------------------------------------

!  Zeroing

     COL_TD = zero ; REFLEC = zero ; SURFACE_TERM = zero
 
!  Zero total reflected contribution (R2_PARTIC) before calculation
!  For Lambertian reflectance, all streams are the same
!  For BRDF, code added 4 May 2009 by R. Spurr

!  Rob Fix 2p4a, 5/26/17. Be careful with energy conservation

      MFAC = HALF * SURFACE_FACTOR

      IF ( DO_INCLUDE_SURFACE ) THEN
        IF ( DO_BRDF_SURFACE ) THEN
          REFLEC = MFAC * BRDF_F(FOURIER)
        ELSE
          REFLEC = MFAC * ALBEDO
        ENDIF
      ENDIF

!  Surface Term------------
!    lowest (surface) boundary with albedo (diffuse radiation terms only)
!    with non-zero albedo, include integrated downward reflectances
!    no albedo, similar code excluding integrated reflectance
!    Add direct beam solution (only to final level)
!    Add thermal emission of ground surface (only to final level)
!mick fix 2/14/2012 - changed treatment of emissivity to be consistent
!                     with LIDORT & VLIDORT

      IF ( DO_INCLUDE_SURFACE ) THEN
        IF ( DO_INCLUDE_DIRECTBEAM ) THEN
          SURFACE_TERM = SURFACE_TERM + DIRECT_BEAM(IPARTIC)
        ENDIF
      ENDIF
      IF ( DO_INCLUDE_SURFEMISS ) THEN
        SURFACE_TERM = SURFACE_TERM + SURFBB * EMISS
      ENDIF

!  --set up Column for solution vector (the "B" as in AX=B)
!  --------------------------------------------------------

! Equivalence: Set up E for tridiagonal system A_l*Y_l-1+B_l*Y_l+D_l*Y_l+1 = E_l

!  Cpt(n) = WUPPER(n,2) ; Cmt(n) = WUPPER(n,1)
!  Cpb(n) = WLOWER(n,2) ; Cmb(n) = WLOWER(n,1)

!  TOA condition
 
!      E(1) = -Cmt(1)
       COL_TD(1) = - WUPPER(1,1)

!  Intermediate layers condition

      do n = 1, nlayers-1
         j1 = 2*n+1 ; j2 = j1 - 1
!         E(j1) = e3(n)*(Cpt(n+1)-Cpb(n))-e1(n)*(Cmt(n+1)-Cmb(n))
!         E(j2) = e2(n+1)*(Cpt(n+1)-Cpb(n))-e4(n+1)*(Cmt(n+1)-Cmb(n)) ! corrected equation
         COL_TD(j1) = emat_td(n,3) * ( WUPPER(2,n+1) - WLOWER(2,n) ) - &
                      emat_td(n,1) * ( WUPPER(1,n+1) - WLOWER(1,n) )
         IF ( CONSSCAT(N+1) ) THEN
!         IF ( CONSSCAT(N) ) THEN
           COL_TD(j2) = WUPPER(2,n+1) - WLOWER(2,n)
         ELSE
           COL_TD(j2) = emat_td(n+1,2) * ( WUPPER(2,n+1) - WLOWER(2,n) ) - &
                        emat_td(n+1,4) * ( WUPPER(1,n+1) - WLOWER(1,n) )
         ENDIF
      enddo

!  Surface condition

      j = ntotal
!      E(j) = S-Cpb(nlayers)+R*Cmb(nlayers)
      COL_TD(j) = SURFACE_TERM - WLOWER(2,nlayers) + REFLEC * WLOWER(1,nlayers)

!      write(*,*)'E(j) at surface = ',SURFACE_TERM,WLOWER(2,nlayers) ,REFLEC,WLOWER(1,nlayers)

!write(0,*) COL_TD(1), SURFACE_TERM, - WLOWER(2,nlayers), REFLEC * WLOWER(1,nlayers), &
!COL_TD(1)+COL_TD(ntotal)

!write(*,*)'check # 18'
!do n = 1, nlayers
!   write(*,*)n,WUPPER(n,1:2),wLOWER(n,1:2)
!enddo

!  Solve tridiagonal system
!  ------------------------

      IF ( NLAYERS .GT. 1 ) THEN

        j = ntotal
        AS(j) = AMAT_TD(j) / BMAT_TD(j)
        DS(j) = COL_TD(j)  / BMAT_TD(j)

        do n = ntotal-1, 1, -1
           X(n) = one/(BMAT_TD(n)-DMAT_TD(n)*AS(n+1))
           AS(n) = AMAT_TD(n)*X(n)
           DS(n) = (COL_TD(n)-DMAT_TD(n)*DS(n+1))*X(n)
        enddo

        Y(1) = DS(1)
        do n = 2, ntotal
           Y(n) = DS(n)-AS(n)*Y(n-1)
        enddo

      ELSE IF ( NLAYERS .EQ. 1 ) THEN

        Y(1) = SELM(1,1) * COL_TD(1) + SELM(1,2) * COL_TD(2)
        Y(2) = SELM(2,1) * COL_TD(1) + SELM(2,2) * COL_TD(2)

      ENDIF

      do n = 1, nlayers
         j1 = 2*n-1 ; j2 = j1 + 1
         Y1(n) = Y(j1)
         Y2(n) = Y(j2)
      enddo

!  Map back to our notation

      do n = 1, nlayers
        if ( consscat(n) ) then
          LCON(N) = Y2(n)
          MCON(N) = Y1(n)  
        endif
      enddo

!write(0,*) LCON(1), MCON(1)
!pause

!  Alternative  Associated quantities
! open(44,file='GammaCheat.dat',status='old')
! do n = 1, nlayers
!    read(44,*)j,mfac,gammad(n),mfac
!enddo
!close(44)

      DO N = 1, NLAYERS
        if ( .not. consscat(n) ) then
!          GAMMA = XPOS(2,N) / XPOS(1,N)
!          GAMMA = GAMMAd(n)
          K1 =   Y1(n) + Y2(n)
          K2 =   Y1(n) - Y2(n)
          LCON_XVEC(1,N) = K2            ; LCON_XVEC(2,N) = K2 * GAMMA(N)
          MCON_XVEC(1,N) = K1 * GAMMA(N) ; MCON_XVEC(2,N) = K1
        endif
      ENDDO

!  Finish

      RETURN
END SUBROUTINE TWOSTREAM_BVP_SOLUTION_TRIDIAG

end module twostream_bvproblem_m
