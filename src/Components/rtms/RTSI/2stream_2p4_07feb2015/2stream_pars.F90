!  P.Castellanos April 2025
!  This is 2STREAM.PARS.


!  File of constants for 2stream model.
!  adapted from VLIDORT version - PC

      MODULE twostream_pars

      IMPLICIT NONE

!  Basic dimensions
!  ================

!  Computational dimensioning
!  --------------------------

!  External thread parameter

      INTEGER, PARAMETER :: MAXTHREADS = 6


!  Maximum number of computational layers

      INTEGER, PARAMETER :: MAXLAYERS = 72
      INTEGER, PARAMETER :: MAXTOTAL = 2*MAXLAYERS

!  Exception handling, maximum number of messages

   INTEGER, PARAMETER :: MAXMESSAGES = 25


!  Geometrical and output parameters
!  ---------------------------------

!  Maximum number of solar zenith angles

      INTEGER, PARAMETER :: MAX_SZANGLES = 2
			    
!  maximum number of user-defined viewing zenith angles

      INTEGER, PARAMETER :: MAX_USER_ANGLES = 2

!  maximum number of user-defined output relative azimuth angles

      INTEGER, PARAMETER :: MAX_USER_RELAZMS = 2

!  Maximum number of Observational Geometries

      INTEGER, PARAMETER :: MAX_USER_OBSGEOMS = 2 !30


!  Surface BRDF dimensioning
!  -------------------------

!  Maximum number of BRDF kernels

      INTEGER, PARAMETER :: MAX_BRDF_KERNELS = 3

!  Maximum number of BRDF parameters per kernel

      INTEGER, PARAMETER :: MAX_BRDF_PARAMETERS = 3

!  Maximum number of azimuth-quadrature streams for BRDF Fourier.
!    5/5/20. Version 2.8.1 Upgrade. MUST BE AN EVEN NUMBER

!     INTEGER, PARAMETER :: MAXSTREAMS_BRDF = 2
      INTEGER, PARAMETER :: MAXSTREAMS_BRDF = 100    ! best


!  Weighting functions
!  -------------------

!  Maximum number of profile/column weighting functions

      INTEGER, PARAMETER :: MAX_ATMOSWFS = 3

!  Maximum number of surface property weighting functions

      INTEGER, PARAMETER :: MAX_SURFACEWFS = 6

!  Maximum number of surface-leaving weighting functions

      INTEGER, PARAMETER :: MAX_SLEAVE_WFS = 1


!  Derived dimensions
!  ==================

!  Copy Beam dimensioning

      INTEGER, PARAMETER :: MAXBEAMS = MAX_SZANGLES

!  Copy viewing zenith angle dimension

      INTEGER, PARAMETER :: MAX_USER_STREAMS = MAX_USER_ANGLES

!  Maximum possible geometries

      INTEGER, PARAMETER :: MAX_GEOMETRIES = &
                            MAX_USER_ANGLES*MAX_USER_RELAZMS*MAX_SZANGLES

!  Numbers
!  ==================
      INTEGER, PARAMETER :: dp = KIND( 1.0D0 )
      INTEGER, parameter :: one = 1.0d0
      REAL(kind=dp), parameter  :: TAYLOR_SMALL = 1.0d-03
!  End of file.

      END MODULE twostream_pars
