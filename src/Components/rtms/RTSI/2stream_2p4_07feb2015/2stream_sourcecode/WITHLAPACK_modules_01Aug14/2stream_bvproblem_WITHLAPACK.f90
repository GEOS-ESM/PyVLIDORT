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

! ###############################################################
! #                                                             #
! # Regular BVP: Subroutines in this Module                     #
! #                                                             #
! #            TWOSTREAM_BVP_MATSETUP_LAPACK                    #
! #            TWOSTREAM_BVP_SOLUTION_LAPACK                    #
! #                                                             #
! #            TWOSTREAM_BVP_MATSETUP_PENTADIAG                 #
! #            TWOSTREAM_BVP_SOLUTION_PENTADIAG                 #
! #                                                             #
! ###############################################################

!  LAPACK routines were reintroduced, 08 April 2014

module twostream_bvproblem_WITHLAPACK_m

   use twostream_LAPACK_m, only : TWOS_DGBTRF, TWOS_DGBTRS 

PUBLIC

contains

SUBROUTINE TWOSTREAM_BVP_MATSETUP_LAPACK &
         ( MAXLAYERS, MAXTOTAL, DO_INCLUDE_SURFACE, FF,             & ! input
           FOURIER, NLAYERS, NTOTAL,                                & ! input
           DO_BRDF_SURFACE, SURFACE_FACTOR, ALBEDO, BRDF_F,         & ! input
           XPOS, XNEG, T_DELT_EIGEN, STREAM_VALUE,                  & ! input
           H_HOMP, H_HOMM, BMAT_ROWMASK, BANDMAT2, IPIVOT, SMAT2,   & ! output
           STATUS, MESSAGE )                                          ! output

      implicit none

!  precision

      INTEGER, PARAMETER :: dp     = KIND( 1.0D0 )

!  input
!  -----

!  Dimensions

      INTEGER, INTENT(IN)        :: MAXLAYERS, MAXTOTAL
      REAL(kind=dp), INTENT(IN)  :: FF

!  control


      LOGICAL, INTENT(IN)        :: DO_INCLUDE_SURFACE
      INTEGER, INTENT(IN)        :: FOURIER
      INTEGER, INTENT(IN)        :: NLAYERS, NTOTAL

!  Surface control

      LOGICAL      , INTENT(IN)  :: DO_BRDF_SURFACE
      REAL(kind=dp), INTENT(IN)  :: SURFACE_FACTOR
      REAL(kind=dp), INTENT(IN)  :: ALBEDO
      REAL(kind=dp), INTENT(IN)  :: BRDF_F(0:1)

!  Eigenvector solutions

      REAL(kind=dp), INTENT(IN)  :: XPOS(2,MAXLAYERS)
      REAL(kind=dp), INTENT(IN)  :: XNEG(2,MAXLAYERS)

!  transmittance factors for +/- eigenvalues

      REAL(kind=dp), INTENT(IN)  :: T_DELT_EIGEN(MAXLAYERS)

!  Stream

      REAL(kind=dp), INTENT(IN)  :: STREAM_VALUE

!  Output
!  ------

!  Downwelling BOA solutions

      REAL(kind=dp), INTENT(OUT) :: H_HOMP
      REAL(kind=dp), INTENT(OUT) :: H_HOMM

!  Initialization of BVP matrix

      INTEGER, INTENT(INOUT)     :: BMAT_ROWMASK(MAXTOTAL,MAXTOTAL)

!  Matrix, Band-matrix for solving BCs

      REAL(kind=dp), INTENT(OUT) :: BANDMAT2(7,MAXTOTAL)

!  square matrix for the single layer case

      REAL(kind=dp), INTENT(OUT) :: SMAT2   (2,2)

!  Pivot matrix

      INTEGER, INTENT(OUT)       :: IPIVOT  (MAXTOTAL)

!  Pentadiagonal matrices
!      REAL(kind=dp), INTENT(OUT) :: MAT(MAXTOTAL,5)
!      REAL(kind=dp), INTENT(OUT) :: ELM (MAXTOTAL,4)
!      REAL(kind=dp), INTENT(OUT) :: SELM (2,2)

!  status

      INTEGER, INTENT(OUT)        :: STATUS
      CHARACTER*(*), INTENT(OUT)  :: MESSAGE

!  local variables
!  ---------------

!  Help

      INTEGER           :: I, N, N1, J, INFO
      INTEGER           :: CP, CM, CEP, CEM, CEP1, CEM1, C0
      REAL(kind=dp)     :: XPNET, XMNET, FACTOR, R2_HOMP, R2_HOMM
      REAL(kind=dp)     :: A, B, C, D, DET
      INTEGER           :: NMIN(MAXTOTAL),NMAX(MAXTOTAL)
      CHARACTER(LEN=3)  :: CI

!  Status

      status = 0
      message = ' '

!  initialize 

      BANDMAT2 = 0.0_dp ; IPIVOT = 0 ; SMAT2 = 0.0_dp

!  Additional setups for the lowest layer
!  For Lambertian reflectance, all streams are the same
!  For BRDF, code added 4 May 2009 by R. Spurr

      H_HOMP = XPOS(1,NLAYERS) * STREAM_VALUE
      H_HOMM = XNEG(1,NLAYERS) * STREAM_VALUE

      R2_HOMP = 0.0d0
      R2_HOMM = 0.0d0
      IF ( DO_INCLUDE_SURFACE ) THEN
        IF ( DO_BRDF_SURFACE ) THEN
          FACTOR = SURFACE_FACTOR * BRDF_F(FOURIER)
        ELSE
          FACTOR = SURFACE_FACTOR * ALBEDO
        ENDIF
        R2_HOMP = FACTOR * H_HOMP
        R2_HOMM = FACTOR * H_HOMM
      ENDIF

!  Inclusion of surface contribution in BV Problem matrix

      IF ( DO_INCLUDE_SURFACE ) THEN
        XPNET = XPOS(2,NLAYERS) - R2_HOMP
        XMNET = XNEG(2,NLAYERS) - R2_HOMM
      ELSE
        XPNET = XPOS(2,NLAYERS)
        XMNET = XNEG(2,NLAYERS)
      ENDIF

!  Initialize Compression Algorithm
!  --------------------------------

      IF ( FOURIER .EQ. 0 ) THEN
        NMIN(1:3) = 1
        DO J = 4, NTOTAL
          NMIN(J) = J - 2
        ENDDO
        DO J = 1, NTOTAL - 2
          NMAX(J) = J + 2
        ENDDO
        NMAX(NTOTAL-1:NTOTAL) = NTOTAL
        DO I = 1, NTOTAL
          DO J = 1, NTOTAL
            IF ( (I.GE.NMIN(J)) .AND. (I.LE.NMAX(J)) ) THEN
              BMAT_ROWMASK(I,J) = 5 + I - J
            ENDIF
          ENDDO
        ENDDO
      ENDIF

!  set up BVP matrix
!  -----------------

!  If Nlayers = 1, special case, ordinary 2x2

      IF ( NLAYERS .EQ. 1 ) THEN

!  Only Top and bottom BC (with surface reflection)

        SMAT2(1,1) = XPOS(1,NLAYERS)
        SMAT2(1,2) = XNEG(1,NLAYERS) * T_DELT_EIGEN(NLAYERS)
        SMAT2(2,1) = XPNET * T_DELT_EIGEN(NLAYERS)
        SMAT2(2,2) = XMNET

!  Invert 2x2 matrix

        A = SMAT2(1,1) ; B = SMAT2(1,2)
        C = SMAT2(2,1) ; D = SMAT2(2,2)
        DET = (A*D - B*C)
        IF ( ABS(DET) .LT. 1.0d-15 ) THEN
          MESSAGE = ' Zero determinant (nlayers=1); '
          STATUS = 1 ; RETURN
        ELSE
          DET = 1.0d0 / DET
          SMAT2(1,1) =   D * DET  ; SMAT2(1,2) = - B * DET 
          SMAT2(2,1) = - C * DET  ; SMAT2(2,2) =   A * DET
        ENDIF

!  General case NLAYERS > 1

      ELSE

!  top BC for layer 1: no downward diffuse radiation

        N = 1
        BANDMAT2(BMAT_ROWMASK(1,1),1)  = XPOS(1,N)
        BANDMAT2(BMAT_ROWMASK(1,2),2)  = XNEG(1,N)*T_DELT_EIGEN(N)

!  intermediate layer boundaries (will not be done if NLAYERS = 1 )

        C0 = - 1
        DO N = 2, NLAYERS
          N1 = N - 1
          C0   = C0 + 2
          DO I = 1, 2
            CM = C0 + I
            CEP  = C0      ; CEM  = CEP + 1
            CEP1 = CEP + 2 ; CEM1 = CEM + 2
            BANDMAT2(BMAT_ROWMASK(CM,CEP),CEP)   =   XPOS(I,N1)*T_DELT_EIGEN(N1)
            BANDMAT2(BMAT_ROWMASK(CM,CEM),CEM)   =   XNEG(I,N1)
            BANDMAT2(BMAT_ROWMASK(CM,CEP1),CEP1) = - XPOS(I,N)
            BANDMAT2(BMAT_ROWMASK(CM,CEM1),CEM1) = - XNEG(I,N) *T_DELT_EIGEN(N)
          ENDDO
        ENDDO

!  bottom BC (with albedo additions if flagged)

        N = NLAYERS
        C0  = C0 + 2  ; CP  = C0 + 1
        CEP = C0     ; CEM = CEP + 1
        BANDMAT2(BMAT_ROWMASK(CP,CEP),CEP) = T_DELT_EIGEN(N) * XPNET
        BANDMAT2(BMAT_ROWMASK(CP,CEM),CEM) = XMNET

!  Scaling

        BANDMAT2 = FF * BANDMAT2

!  LAPACK LU-decomposition for band matrix

        CALL TWOS_DGBTRF ( NTOTAL, NTOTAL, 2, 2, BANDMAT2, 7, IPIVOT, INFO )

!  Exception handling

        IF ( INFO .GT. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO
          MESSAGE = 'Singular matrix, u(i,i)=0, for i = '//CI//'; TWOS_DGBTRF call '
          STATUS = 1 ; RETURN
        ELSE IF ( INFO .LT. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO
          MESSAGE = 'argument i illegal value, for i = '//CI//'; TWOS_DGBTRF call '
          STATUS = 1 ; RETURN
        ENDIF

      ENDIF

!  finish

      RETURN
END SUBROUTINE TWOSTREAM_BVP_MATSETUP_LAPACK

!

SUBROUTINE TWOSTREAM_BVP_SOLUTION_LAPACK &
      ( MAXLAYERS, MAXBEAMS, MAXTOTAL, FF,               & ! Dimensions
        DO_INCLUDE_SURFACE, DO_INCLUDE_DIRECTBEAM,       & ! inputs
        DO_INCLUDE_SURFEMISS, DO_BRDF_SURFACE,           & ! inputs
        FOURIER, IPARTIC, NLAYERS, NTOTAL,               & ! inputs
        SURFACE_FACTOR, ALBEDO, BRDF_F, EMISS, SURFBB,   & ! inputs
        DIRECT_BEAM, XPOS, XNEG, WUPPER, WLOWER,         & ! inputs
        STREAM_VALUE, BANDMAT2, SMAT2, IPIVOT,           & ! inputs
        H_PARTIC, LCON, MCON, LCON_XVEC, MCON_XVEC,      & ! output
        STATUS, MESSAGE )                                  ! output

      implicit none

!  Precision

      INTEGER, PARAMETER :: dp     = KIND( 1.0D0 )

!  Input arguments
!  ---------------

!  Dimensions

      INTEGER, INTENT(IN)        :: MAXBEAMS, MAXLAYERS, MAXTOTAL
      REAL(kind=dp), INTENT(IN)  :: FF

!  Inclusion flags

      LOGICAL, INTENT(IN)        :: DO_INCLUDE_DIRECTBEAM
      LOGICAL, INTENT(IN)        :: DO_INCLUDE_SURFACE
      LOGICAL, INTENT(IN)        :: DO_INCLUDE_SURFEMISS

!  Surface control

      LOGICAL      , INTENT(IN)  :: DO_BRDF_SURFACE
      REAL(kind=dp), INTENT(IN)  :: SURFACE_FACTOR
      REAL(kind=dp), INTENT(IN)  :: ALBEDO
      REAL(kind=dp), INTENT(IN)  :: BRDF_F(0:1)
      REAL(kind=dp), INTENT(IN)  :: EMISS, SURFBB

!  Fourier component and beam number

      INTEGER, INTENT(IN)        :: FOURIER, IPARTIC

!  Numbers

      INTEGER, INTENT(IN)        :: NLAYERS, NTOTAL

!  Direct beam

      REAL(kind=dp), INTENT(IN)  :: DIRECT_BEAM ( MAXBEAMS )

!  Eigenvector solutions

      REAL(kind=dp), INTENT(IN)  :: XPOS(2,MAXLAYERS)
      REAL(kind=dp), INTENT(IN)  :: XNEG(2,MAXLAYERS)

!  Particular solutions

      REAL(kind=dp), INTENT(IN)  :: WLOWER ( 2, MAXLAYERS )
      REAL(kind=dp), INTENT(IN)  :: WUPPER ( 2, MAXLAYERS )

!  Stream

      REAL(kind=dp), INTENT(IN)  :: STREAM_VALUE

!  Matrix, Band-matrix for solving BCs

      REAL(kind=dp), INTENT(INOUT) :: BANDMAT2(7,MAXTOTAL)

!  square matrix for the single layer case

      REAL(kind=dp), INTENT(INOUT) :: SMAT2   (2,2)

!  Pivot matrix

      INTEGER, INTENT(INOUT)       :: IPIVOT  (MAXTOTAL)

!  Pentadiagonal Matrix entries for solving BCs, Elimination mattrices
!      REAL(kind=dp), INTENT(IN)  :: MAT(MAXTOTAL,5)
!      REAL(kind=dp), INTENT(IN)  :: ELM (MAXTOTAL,4)
!      REAL(kind=dp), INTENT(IN)  :: SELM (2,2)

!  Output
!  ------

!  Downwelling BOA solution

      REAL(kind=dp), INTENT(OUT) :: H_PARTIC

!  Solution constants of integration, and related quantities

      REAL(kind=dp), INTENT(OUT) :: LCON(MAXLAYERS)
      REAL(kind=dp), INTENT(OUT) :: MCON(MAXLAYERS)

      REAL(kind=dp), INTENT(OUT) :: LCON_XVEC(2,MAXLAYERS)
      REAL(kind=dp), INTENT(OUT) :: MCON_XVEC(2,MAXLAYERS)

!  status

      INTEGER, INTENT(OUT)        :: STATUS
      CHARACTER*(*), INTENT(OUT)  :: MESSAGE

!  Local variables
!  ---------------

!  Column vectors for solving BCs. Not saved.

      REAL(kind=dp)       :: COL2   (MAXTOTAL,1)
      REAL(kind=dp)       :: SCOL2  (2)

!  Other variables

      INTEGER             :: N, N1, I, NM, NP, INFO
      REAL(kind=dp)       :: FACTOR, R2_PARTIC, A, B
      character*3         :: CI

!  Status

      status = 0
      message = ' '

!  --Additional setups for the bottom layer
!  ----------------------------------------

!  Zero total reflected contributions
!  For Lambertian reflectance, all streams are the same
!  For BRDF, code added 4 May 2009 by R. Spurr

      H_PARTIC = WLOWER(1,NLAYERS) * STREAM_VALUE
      R2_PARTIC = 0.0d0
      IF ( DO_INCLUDE_SURFACE ) THEN
        IF ( DO_BRDF_SURFACE ) THEN
          FACTOR = SURFACE_FACTOR * BRDF_F(FOURIER)
        ELSE
          FACTOR = SURFACE_FACTOR * ALBEDO
        ENDIF
        R2_PARTIC = H_PARTIC * FACTOR
      ENDIF

!  --set up Column for solution vector (the "B" as in AX=B)
!  --------------------------------------------------------

!  Pentadiagonal case

      IF ( NLAYERS .GT. 1 ) THEN

!  zero column vector

        COL2 = 0.0d0

!  Upper boundary for layer 1: no downward diffuse radiation

        COL2(1,1)   = - WUPPER(1,1)

!  intermediate layer boundaries

        DO N = 2, NLAYERS
          N1 = N - 1
          NM = 2*N1 ; NP = NM + 1
          COL2(NM,1) = WUPPER(1,N) - WLOWER(1,N1)
          COL2(NP,1) = WUPPER(2,N) - WLOWER(2,N1)
        ENDDO

!  lowest (surface) boundary with albedo (diffuse + direct)

        COL2(NTOTAL,1) = - WLOWER(2,NLAYERS)
        IF ( DO_INCLUDE_SURFACE ) THEN
          COL2(NTOTAL,1) = COL2(NTOTAL,1) + R2_PARTIC
          IF ( DO_INCLUDE_DIRECTBEAM ) THEN
            COL2(NTOTAL,1) = COL2(NTOTAL,1) + DIRECT_BEAM(IPARTIC)
          ENDIF
        ENDIF

!  Add thermal emission of ground surface (only to final level)
!mick fix 2/14/2012 - changed treatment of emissivity to be consistent with LIDORT & VLIDORT

        !IF ( DO_INCLUDE_SURFEMISS ) THEN
        !  IF ( DO_BRDF_SURFACE ) THEN
        !    LOCAL_EMISS = EMISS
        !  ELSE
        !    LOCAL_EMISS = 1.0_dp - ALBEDO
        !  ENDIF
        !  COL(NTOTAL) = COL(NTOTAL) + SURFBB * local_emiss
        !ENDIF

        IF ( DO_INCLUDE_SURFEMISS ) THEN
          COL2(NTOTAL,1) = COL2(NTOTAL,1) + SURFBB * EMISS
        ENDIF

!  Scaling

        COL2 = FF * COL2

!  LAPACK substitution (DGBTRS) using RHS column vector COL2

        CALL TWOS_DGBTRS( 'n', NTOTAL, 2, 2, 1, BANDMAT2, 7, IPIVOT, &
                     COL2, MAXTOTAL, INFO )

!  Exception handling

        IF ( INFO .LT. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO
          MESSAGE = 'argument i illegal value, for i = '//CI//'TWOS_DGBTRS call'
          STATUS = 1 ; RETURN
         ENDIF

!  Set integration constants LCON and MCON for -/+ eigensolutions, all layers

        DO N = 1, NLAYERS
          NM = 2*N-1 ; NP = NM + 1
          LCON(N) = COL2(NM,1)
          MCON(N) = COL2(NP,1)
        ENDDO

!  If Nlayers = 1, special case

      ELSE IF ( NLAYERS .EQ. 1 ) THEN

!  Upper boundary for layer 1: no downward diffuse radiation
!  lowest (surface) boundary with albedo (diffuse radiation terms only)
!  with non-zero albedo, include integrated downward reflectances
!  no albedo, similar code excluding integrated reflectance
!  Add direct beam solution (only to final level)
!  Add thermal emission of ground surface (only to final level)

        SCOL2 = 0.0d0
        SCOL2(1:2) = - WUPPER(1:2,1)
        IF ( DO_INCLUDE_SURFACE ) THEN
          SCOL2(2) = SCOL2(2)  + R2_PARTIC
          IF ( DO_INCLUDE_DIRECTBEAM ) THEN
            SCOL2(2) = SCOL2(2) + DIRECT_BEAM(IPARTIC)
          ENDIF
        ENDIF
        IF ( DO_INCLUDE_SURFEMISS ) THEN
          SCOL2(2) = SCOL2(2) + SURFBB * EMISS
        ENDIF

!  Solve the boundary problem: No compression, Single Layer only

        A = SCOL2(1) ; B = SCOL2(2)
        SCOL2(1) = SMAT2(1,1) * A + smat2(1,2) * B
        SCOL2(2) = SMAT2(2,1) * A + smat2(2,2) * B
        LCON(1) = SCOL2(1)
        MCON(1) = SCOL2(2)

      ENDIF

!  Associated quantities
!  ---------------------

      DO N = 1, NLAYERS
        DO I = 1, 2
          LCON_XVEC(I,N) = LCON(N)*XPOS(I,N)
          MCON_XVEC(I,N) = MCON(N)*XNEG(I,N)
        ENDDO
      ENDDO

!  Finish

      RETURN
END SUBROUTINE TWOSTREAM_BVP_SOLUTION_LAPACK



SUBROUTINE TWOSTREAM_BVP_MATSETUP_PENTADIAG &
         ( MAXLAYERS, MAXTOTAL, FF, DO_INVERSE,                     & ! Dimensions
           DO_INCLUDE_SURFACE, FOURIER_COMPONENT, NLAYERS, NTOTAL,  & ! input
           DO_BRDF_SURFACE, SURFACE_FACTOR, ALBEDO, BRDF_F,         & ! input
           XPOS, XNEG, T_DELT_EIGEN, STREAM_VALUE,                  & ! input
           H_HOMP, H_HOMM, MAT, ELM, SELM,                          & ! output
           STATUS, MESSAGE )                                          ! output

      implicit none

!  precision

      INTEGER, PARAMETER :: dp     = KIND( 1.0D0 )

!  input
!  -----

!  Dimensions

      INTEGER, INTENT(IN)        :: MAXLAYERS, MAXTOTAL

!  Inverse control

      LOGICAL      , INTENT(IN)  :: DO_INVERSE
      REAL(kind=dp), INTENT(IN)  :: FF

!  control

      LOGICAL, INTENT(IN)        :: DO_INCLUDE_SURFACE
      INTEGER, INTENT(IN)        :: FOURIER_COMPONENT
      INTEGER, INTENT(IN)        :: NLAYERS, NTOTAL

!  Surface control

      LOGICAL      , INTENT(IN)  :: DO_BRDF_SURFACE
      REAL(kind=dp), INTENT(IN)  :: SURFACE_FACTOR
      REAL(kind=dp), INTENT(IN)  :: ALBEDO
      REAL(kind=dp), INTENT(IN)  :: BRDF_F(0:1)

!  Eigenvector solutions

      REAL(kind=dp), INTENT(IN)  :: XPOS(2,MAXLAYERS)
      REAL(kind=dp), INTENT(IN)  :: XNEG(2,MAXLAYERS)

!  transmittance factors for +/- eigenvalues

      REAL(kind=dp), INTENT(IN)  :: T_DELT_EIGEN(MAXLAYERS)

!  Stream

      REAL(kind=dp), INTENT(IN)  :: STREAM_VALUE

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
      REAL(kind=dp)     :: XPNET, XMNET, FACTOR, BET, DEN, R2_HOMP, R2_HOMM
!      REAL(kind=dp)     :: IMAT(MAXTOTAL,5)
      CHARACTER(LEN=3)  :: CI

!  Stability check value

      REAL(kind=dp)     :: SMALLNUM=1.0D-20

!  Status

      status = 0
      message = ' '

!  Additional setups for the lowest layer
!  For Lambertian reflectance, all streams are the same
!  For BRDF, code added 4 May 2009 by R. Spurr

      H_HOMP = XPOS(1,NLAYERS) * STREAM_VALUE
      H_HOMM = XNEG(1,NLAYERS) * STREAM_VALUE

      R2_HOMP = 0.0d0
      R2_HOMM = 0.0d0
      IF ( DO_INCLUDE_SURFACE ) THEN
        IF ( DO_BRDF_SURFACE ) THEN
          FACTOR = SURFACE_FACTOR * BRDF_F(FOURIER_COMPONENT)
        ELSE
          FACTOR = SURFACE_FACTOR * ALBEDO
        ENDIF
        R2_HOMP = FACTOR * H_HOMP
        R2_HOMM = FACTOR * H_HOMM
      ENDIF

!  Inclusion of surface contribution in BV Problem matrix

      IF ( DO_INCLUDE_SURFACE ) THEN
        XPNET = XPOS(2,NLAYERS) - R2_HOMP
        XMNET = XNEG(2,NLAYERS) - R2_HOMM
      ELSE
        XPNET = XPOS(2,NLAYERS)
        XMNET = XNEG(2,NLAYERS)
      ENDIF

!  set up BVP matrix
!  -----------------

!  If Nlayers = 1, special case, ordinary 2x2

      IF ( NLAYERS .EQ. 1 ) THEN

!  Only Top and bottom BC (with surface reflection)

        SMAT(1,1) = XPOS(1,NLAYERS)
        SMAT(1,2) = XNEG(1,NLAYERS) * T_DELT_EIGEN(NLAYERS)
        SMAT(2,1) = XPNET * T_DELT_EIGEN(NLAYERS)
        SMAT(2,2) = XMNET

!  If NLAYERS > 1, set up Pentadiagonal matrix

      ELSE

!  Zero for both Fourier components (Important bug!)

        MAT = 0.0d0
        ELM = 0.0d0

!  top BC for layer 1: no downward diffuse radiation
!  intermediate layer boundaries
!  bottom BC (including surface reflected term)

!  Regular set

        if ( .NOT. do_inverse ) then
          MAT(1,3)  = XPOS(1,1)
          MAT(1,4)  = XNEG(1,1) * T_DELT_EIGEN(1)
          DO N = 2, NLAYERS
            N1 =  N - 1
            NM = 2*N1
            NP = NM + 1
            MAT(NM,2) =   XPOS(1,N1) * T_DELT_EIGEN(N1)
            MAT(NM,3) =   XNEG(1,N1)
            MAT(NM,4) = - XPOS(1,N)
            MAT(NM,5) = - XNEG(1,N)  * T_DELT_EIGEN(N)
            MAT(NP,1) =   XPOS(2,N1) * T_DELT_EIGEN(N1)
            MAT(NP,2) =   XNEG(2,N1)
            MAT(NP,3) = - XPOS(2,N)
            MAT(NP,4) = - XNEG(2,N)  * T_DELT_EIGEN(N)
          ENDDO
          MAT(NTOTAL,2) = XPNET * T_DELT_EIGEN(NLAYERS)
          MAT(NTOTAL,3) = XMNET
        endif

!  Inverted set

        if ( do_inverse ) then
          MAT(1,3)  = XMNET
          MAT(1,4)  = XPNET * T_DELT_EIGEN(NLAYERS)
          DO N = 2, NLAYERS
            N1 =  N - 1
            INM = NTOTAL - 2*N1  ; INP = INM + 1
            MAT(INM,2) =  - XNEG(2,N)  * T_DELT_EIGEN(N)
            MAT(INM,3) =  - XPOS(2,N)
            MAT(INM,4) =    XNEG(2,N1)
            MAT(INM,5) =    XPOS(2,N1) * T_DELT_EIGEN(N1)
            MAT(INP,1) =  - XNEG(1,N)  * T_DELT_EIGEN(N)
            MAT(INP,2) =  - XPOS(1,N)
            MAT(INP,3) =    XNEG(1,N1)
            MAT(INP,4) =    XPOS(1,N1) * T_DELT_EIGEN(N1)
          ENDDO
          MAT(NTOTAL,2) = XNEG(1,1) * T_DELT_EIGEN(1)
          MAT(NTOTAL,3) = XPOS(1,1)
        endif

!  original code

!        if ( do_inverse ) then
!          IMAT = 0.0d0
!          IMAT(1,3)  = MAT(NTOTAL,3)
!          IMAT(1,4)  = MAT(NTOTAL,2) 
!          DO N = 2, NLAYERS
!            NM  = 2*N - 2          ; NP  = NM + 1
!            INP = NTOTAL + 1 - NM  ; INM = INP - 1
!            IMAT(INM,2) = MAT(NP,4) 
!            IMAT(INM,3) = MAT(NP,3)
!            IMAT(INM,4) = MAT(NP,2)
!            IMAT(INM,5) = MAT(NP,1)
!            IMAT(INP,1) = MAT(NM,5)
!            IMAT(INP,2) = MAT(NM,4)
!            IMAT(INP,3) = MAT(NM,3)
!            IMAT(INP,4) = MAT(NM,2)
!          ENDDO
!          IMAT(NTOTAL,2) = MAT(1,4)
!          IMAT(NTOTAL,3) = MAT(1,3)
!          MAT = IMAT
!        endif

      ENDIF

!  Scaling

      MAT = FF * MAT

!  Elimination of BVP pentadiagonal matrix
!  ---------------------------------------

      IF ( NLAYERS .GT. 1 ) THEN

!  Row 1

        ELM(1,4) = 0.0d0
        ELM(1,3) = 1.0d0/MAT(1,3)
        ELM(1,1) = - MAT(1,4) * ELM(1,3)
        ELM(1,2) = - MAT(1,5) * ELM(1,3)

!  Row 2; includes first check for singularity

        ELM(2,4) = 0.0d0
        bet = MAT(2,3) + MAT(2,2)*ELM(1,1)
        IF ( DABS(Bet) .LT. SMALLNUM ) THEN 
          message = 'Singularity in Pentadiagonal Matrix, Row #  2'
          status = 1
          return
        endif
        bet = -1.0d0/bet
        ELM(2,1) = (MAT(2,4) + MAT(2,2)*ELM(1,2)) * bet
        ELM(2,2) = MAT(2,5) * bet
        ELM(2,3) = bet
        ELM(2,4) = 0.0d0

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
          den = - 1.0d0 / den
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
        DEN = 1.0d0 / DEN
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
      ( MAXLAYERS, MAXBEAMS, MAXTOTAL, FF, DO_INVERSE,   & ! Dimensions
        DO_INCLUDE_SURFACE, DO_INCLUDE_DIRECTBEAM,       & ! inputs
        DO_INCLUDE_SURFEMISS, DO_BRDF_SURFACE,           & ! inputs
        FOURIER, IPARTIC, NLAYERS, NTOTAL,               & ! inputs
        SURFACE_FACTOR, ALBEDO, BRDF_F, EMISS, SURFBB,   & ! inputs
        DIRECT_BEAM, XPOS, XNEG, WUPPER, WLOWER,         & ! inputs
        STREAM_VALUE, MAT, ELM, SELM,                    & ! inputs
        H_PARTIC, LCON, MCON, LCON_XVEC, MCON_XVEC  )      ! Output

      implicit none

!  Precision

      INTEGER, PARAMETER :: dp     = KIND( 1.0D0 )

!  Input arguments
!  ---------------

!  Dimensions

      INTEGER, INTENT(IN)        :: MAXBEAMS, MAXLAYERS, MAXTOTAL

!  Inverse control

      LOGICAL      , INTENT(IN)  :: DO_INVERSE
      REAL(kind=dp), INTENT(IN)  :: FF

!  Inclusion flags

      LOGICAL, INTENT(IN)        :: DO_INCLUDE_DIRECTBEAM
      LOGICAL, INTENT(IN)        :: DO_INCLUDE_SURFACE
      LOGICAL, INTENT(IN)        :: DO_INCLUDE_SURFEMISS

!  Surface control

      LOGICAL      , INTENT(IN)  :: DO_BRDF_SURFACE
      REAL(kind=dp), INTENT(IN)  :: SURFACE_FACTOR
      REAL(kind=dp), INTENT(IN)  :: ALBEDO
      REAL(kind=dp), INTENT(IN)  :: BRDF_F(0:1)
      REAL(kind=dp), INTENT(IN)  :: EMISS, SURFBB

!  Fourier component and beam number

      INTEGER, INTENT(IN)        :: FOURIER, IPARTIC

!  Numbers

      INTEGER, INTENT(IN)        :: NLAYERS, NTOTAL

!  Direct beam

      REAL(kind=dp), INTENT(IN)  :: DIRECT_BEAM ( MAXBEAMS )

!  Eigenvector solutions

      REAL(kind=dp), INTENT(IN)  :: XPOS(2,MAXLAYERS)
      REAL(kind=dp), INTENT(IN)  :: XNEG(2,MAXLAYERS)

!  Particular solutions

      REAL(kind=dp), INTENT(IN)  :: WLOWER ( 2, MAXLAYERS )
      REAL(kind=dp), INTENT(IN)  :: WUPPER ( 2, MAXLAYERS )

!  Stream

      REAL(kind=dp), INTENT(IN)  :: STREAM_VALUE

!  Pentadiagonal Matrix entries for solving BCs

      REAL(kind=dp), INTENT(IN)  :: MAT(MAXTOTAL,5)

!  Pentadiagonal elimination matrix

      REAL(kind=dp), INTENT(IN)  :: ELM (MAXTOTAL,4)

!  Single layer elimination matrix

      REAL(kind=dp), INTENT(IN)  :: SELM (2,2)

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

      INTEGER             :: N, N1, I, NM, NP, INP, INM, NI
      REAL(kind=dp)       :: FACTOR, DEN, R2_PARTIC, TOA_TERM, SURFACE_TERM
      REAL(kind=dp)       :: NEW_SCOL1  !, ICOL(MAXTOTAL)

!  --Additional setups for the bottom layer
!  ----------------------------------------

!  Zero total reflected contributions
!  For Lambertian reflectance, all streams are the same
!  For BRDF, code added 4 May 2009 by R. Spurr

      H_PARTIC = WLOWER(1,NLAYERS) * STREAM_VALUE
      R2_PARTIC = 0.0d0
      IF ( DO_INCLUDE_SURFACE ) THEN
        IF ( DO_BRDF_SURFACE ) THEN
          FACTOR = SURFACE_FACTOR * BRDF_F(FOURIER)
        ELSE
          FACTOR = SURFACE_FACTOR * ALBEDO
        ENDIF
        R2_PARTIC = H_PARTIC * FACTOR
      ENDIF

!  --set up Column for solution vector (the "B" as in AX=B)
!  --------------------------------------------------------

!  Pentadiagonal case

      IF ( NLAYERS .GT. 1 ) THEN

!  zero column vector

        COL(1:NTOTAL) = 0.0d0

!  Upper boundary for layer 1: no downward diffuse radiation

        TOA_TERM  = - WUPPER(1,1)

!  Surface Term------------
!    lowest (surface) boundary with albedo (diffuse + direct)
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

!  vector

        if ( do_inverse ) then
          COL(1)  = SURFACE_TERM
          DO N = 2, NLAYERS
            N1 = N - 1
            INM = NTOTAL - 2*N1
            INP = INM + 1
            COL(INP) = WUPPER(1,N) - WLOWER(1,N1)
            COL(INM) = WUPPER(2,N) - WLOWER(2,N1)
          ENDDO
          COL(NTOTAL) = TOA_TERM
        else
          COL(1)  = TOA_TERM
          DO N = 2, NLAYERS
            N1 = N - 1
            NM = 2*N1
            NP = NM + 1
            COL(NM) = WUPPER(1,N) - WLOWER(1,N1)
            COL(NP) = WUPPER(2,N) - WLOWER(2,N1)
          ENDDO
          COL(NTOTAL) = SURFACE_TERM
        endif

!  Scaling

        COL = FF * COL

!  If Nlayers = 1, special case

      ELSE IF ( NLAYERS .EQ. 1 ) THEN

!  Upper boundary for layer 1: no downward diffuse radiation
!  lowest (surface) boundary with albedo (diffuse radiation terms only)
!  with non-zero albedo, include integrated downward reflectances
!  no albedo, similar code excluding integrated reflectance
!  Add direct beam solution (only to final level)

        SCOL(1) = - WUPPER(1,1)
        SCOL(2) = - WLOWER(2,1)
        IF ( DO_INCLUDE_SURFACE ) THEN
          SCOL(2) = SCOL(2)  + R2_PARTIC
          IF ( DO_INCLUDE_DIRECTBEAM ) THEN
            SCOL(2) = SCOL(2) + DIRECT_BEAM(IPARTIC)
          ENDIF
        ENDIF

!  Add thermal emission of ground surface (only to final level)

!mick fix 2/14/2012 - changed treatment of emissivity to be consistent
!                     with LIDORT & VLIDORT

        !IF ( DO_INCLUDE_SURFEMISS ) THEN
        !  IF ( DO_BRDF_SURFACE ) THEN
        !    LOCAL_EMISS = EMISS
        !  ELSE
        !    LOCAL_EMISS = 1.0_dp - ALBEDO
        !  ENDIF
        !  SCOL(2) = SCOL(2) + SURFBB * local_emiss
        !ENDIF

        IF ( DO_INCLUDE_SURFEMISS ) THEN
          SCOL(2) = SCOL(2) + SURFBB * EMISS
        ENDIF

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
           DO N = 1, NLAYERS
              NI = NLAYERS + 1 - N
              INP = 2*NI
              INM = INP - 1
              LCON(N) = COL(INP)
              MCON(N) = COL(INM)
           ENDDO
        else
           DO N = 1, NLAYERS
              NM = 2*N-1
              NP = NM + 1
              LCON(N) = COL(NM)
              MCON(N) = COL(NP)
           ENDDO
        endif

!  Solve the boundary problem: Single Layer only

      ELSE IF ( NLAYERS .EQ. 1 ) THEN
!mick fix 1/9/2011 - defined NEW_SCOL1 so a MODIFIED version of SCOL(1)
!                    is not being used in computing SCOL(2)
        NEW_SCOL1 = SELM(1,1) * SCOL(1) + SELM(1,2) * SCOL(2)
        SCOL(2)   = SELM(2,1) * SCOL(1) + SELM(2,2) * SCOL(2)

        LCON(1) = NEW_SCOL1
        MCON(1) = SCOL(2)
      ENDIF

!  debug
!      if ( fourier.eq.0 .and.ipartic.eq.1) then
!        do n = 1, nlayers
!          write(*,'(i3,1p2e24.12)')n,LCON(N),MCON(N)
!        enddo
!      endif

!  Associated quantities
!  ---------------------

      DO N = 1, NLAYERS
        DO I = 1, 2
          LCON_XVEC(I,N) = LCON(N)*XPOS(I,N)
          MCON_XVEC(I,N) = MCON(N)*XNEG(I,N)
        ENDDO
      ENDDO

!  Finish

      RETURN
END SUBROUTINE TWOSTREAM_BVP_SOLUTION_PENTADIAG

!  End module

end module twostream_bvproblem_WITHLAPACK_m
