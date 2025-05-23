
! ###############################################################
! #                                                             #
! #                       VLIDORT_2p8p3                         #
! #                                                             #
! #  Vectorized LInearized Discrete Ordinate Radiative Transfer #
! #  -          --         -        -        -         -        #
! #                                                             #
! ###############################################################

! ###############################################################
! #                                                             #
! #  Authors :     Robert. J. D. Spurr (1)                      #
! #                Matt Christi                                 #
! #                                                             #
! #  Address (1) : RT Solutions, inc.                           #
! #                9 Channing Street                            #
! #                Cambridge, MA 02138, USA                     #
! #                                                             #
! #  Tel:          (617) 492 1183                               #
! #  Email :       rtsolutions@verizon.net                      #
! #                                                             #
! #  This Version :   VLIDORT_2p8p3                             #
! #  Release Date :   31 March 2021                             #
! #                                                             #
! #  Previous VLIDORT Versions under Standard GPL 3.0:          #
! #  ------------------------------------------------           #
! #                                                             #
! #      2.7   F90, released        August 2014                 #
! #      2.8   F90, released        May    2017                 #
! #      2.8.1 F90, released        August 2019                 # 
! #      2.8.2 F90, limited release May    2020                 # 
! #                                                             #
! #  Features Summary of Recent VLIDORT Versions:               #
! #  -------------------------------------------                #
! #                                                             #
! #      NEW: TOTAL COLUMN JACOBIANS         (2.4)              #
! #      NEW: BPDF Land-surface KERNELS      (2.4R)             #
! #      NEW: Thermal Emission Treatment     (2.4RT)            #
! #      Consolidated BRDF treatment         (2.4RTC)           #
! #      f77/f90 Release                     (2.5)              #
! #      External SS / New I/O Structures    (2.6)              #
! #                                                             #
! #      SURFACE-LEAVING / BRDF-SCALING      (2.7)              #
! #      TAYLOR Series / OMP THREADSAFE      (2.7)              #
! #      New Water-Leaving Treatment         (2.8)              #
! #      LBBF & BRDF-Telescoping, enabled    (2.8)              #
! #      Several Performance Enhancements    (2.8)              #
! #      Water-leaving coupled code          (2.8.1)            #
! #      Planetary problem, media properties (2.8.1)            #
! #      Doublet geometry post-processing    (2.8.2)            #
! #      Reduction zeroing, dynamic memory   (2.8.2)            #
! #                                                             #
! #  Features Summary of This VLIDORT Version                   #
! #  ----------------------------------------                   #
! #                                                             #
! #   2.8.3, released 31 March 2021.                            #
! #     ==> Green's function RT solutions (Nstokes = 1 or 3)    #
! #     ==> Sphericity Corrections using MS source terms        #
! #     ==> BRDF upgrades, including new snow reflectance       #
! #     ==> SLEAVE Upgrades, extended water-leaving treatment   #
! #                                                             #
! ###############################################################

! ###################################################################
! #                                                                 #
! # This is Version 2.8.3 of the VLIDORT_2p8 software library.      #
! # This library comes with the Standard GNU General Public License,#
! # Version 3.0, 29 June 2007. Please read this license carefully.  #
! #                                                                 #
! #      VLIDORT Copyright (c) 2003-2021.                           #
! #          Robert Spurr, RT Solutions, Inc.                       #
! #          9 Channing Street, Cambridge, MA 02138, USA.           #
! #                                                                 #
! # This file is part of VLIDORT_2p8p3 ( Version 2.8.3 )            #
! #                                                                 #
! # VLIDORT_2p8p3 is free software: you can redistribute it         #
! # and/or modify it under the terms of the Standard GNU GPL        #
! # (General Public License) as published by the Free Software      #
! # Foundation, either version 3.0 of the License, or any           #
! # later version.                                                  #
! #                                                                 #
! # VLIDORT_2p8p3 is distributed in the hope that it will be        #
! # useful, but WITHOUT ANY WARRANTY; without even the implied      #
! # warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR         #
! # PURPOSE. See the Standard GNU General Public License (GPL)      #
! # for more details.                                               #
! #                                                                 #
! # You should have received a copy of the Standard GNU General     #
! # Public License (GPL) Version 3.0, along with the VLIDORT_2p8p3  #
! # code package. If not, see <http://www.gnu.org/licenses/>.       #
! #                                                                 #
! ###################################################################

!  This is VLIDORT.PARS.

!  File of constants for VLIDORT model.

      MODULE vlidort_pars_m

      IMPLICIT NONE

!  Real number type definitions

      INTEGER, PARAMETER :: VLIDORT_SPKIND = SELECTED_REAL_KIND(6)
      INTEGER, PARAMETER :: VLIDORT_DPKIND = SELECTED_REAL_KIND(15)
      INTEGER, PARAMETER :: FPK = VLIDORT_DPKIND

!  RTS use

      INTEGER :: SDU,LDU

!  Version number
!  ==============

      CHARACTER (LEN=3), PARAMETER :: VLIDORT_VERSION_NUMBER = '2.8'

!  File i/o unit numbers
!  ======================


      INTEGER, PARAMETER ::       VLIDORT_INUNIT   = 21
      INTEGER, PARAMETER ::       VLIDORT_SCENUNIT = 22
      INTEGER, PARAMETER ::       VLIDORT_FUNIT    = 23
      INTEGER, PARAMETER ::       VLIDORT_RESUNIT  = 24
      INTEGER, PARAMETER ::       VLIDORT_ERRUNIT  = 25

!  Basic dimensions
!  ================

!  Computational dimensioning
!  --------------------------

!  Number of computational streams in the half-space

      INTEGER, PARAMETER :: MAXSTREAMS = 25

!  Maximum number of computational layers

      INTEGER, PARAMETER :: MAXLAYERS = 72

!  Maximum number of fine layers used in single scattering corrections

      INTEGER, PARAMETER :: MAXFINELAYERS = 4

!  Maximum number of input moments.
!    (Use full range for exact single scatter calculations)

      INTEGER, PARAMETER :: MAXMOMENTS_INPUT = 1000

!  Max number of thermal coefficients
!   --------- New for Version 2.4RT -----------------

      INTEGER, PARAMETER :: MAX_THERMAL_COEFFS = 2

!  Geometrical and output parameters
!  ---------------------------------

!  Maximum number of solar zenith angles

!      INTEGER, PARAMETER :: MAX_SZANGLES = 21 !14 !30  
      INTEGER, PARAMETER :: MAX_SZANGLES = 2
			    
!  maximum number of user-defined viewing zenith angles

!      INTEGER, PARAMETER :: MAX_USER_VZANGLES = 21 ! 15 !30
      INTEGER, PARAMETER :: MAX_USER_VZANGLES = 2

!  maximum number of user-defined output relative azimuth angles

!      INTEGER, PARAMETER :: MAX_USER_RELAZMS = 21 !14 !30 
      INTEGER, PARAMETER :: MAX_USER_RELAZMS = 2

!  Maximum number of Observational Geometries
!   New parameter, 25 October 2012

      INTEGER, PARAMETER :: MAX_USER_OBSGEOMS = 2 !30

!  Maximum number of output optical depths

      INTEGER, PARAMETER :: MAX_USER_LEVELS = 2

!  Maximum number of output optical depths away from layer boundaries
!   This must be less than or equal to the previous entry

      INTEGER, PARAMETER :: MAX_PARTLAYERS = 2

!  Version 2p7. Maximum number of Terms for Taylor series expansions
!    If you are retaining contributions of order EPS^n, 
!    then you need at least n+2 Taylor terms

   INTEGER, PARAMETER :: MAX_TAYLOR_TERMS = 7

!  Fixed parameters
!  ----------------

!  Number of Stokes parameters

      INTEGER, PARAMETER :: MAXSTOKES = 4

!  Two directions (Up and Down)

      INTEGER, PARAMETER :: MAX_DIRECTIONS = 2

!  Surface BRDF dimensioning
!  -------------------------

!  Maximum number of BRDF kernels
!    ** Increased number of kernels to 4, Version 2.8

      INTEGER, PARAMETER :: MAX_BRDF_KERNELS = 4

!  Maximum number of BRDF parameters per kernel
!    ** Increased number of parameters to 4, Version 2.8

      INTEGER, PARAMETER :: MAX_BRDF_PARAMETERS = 4

!  Maximum number of azimuth-quadrature streams for BRDF Fourier.
!    5/5/20. Version 2.8.1 Upgrade. MUST BE AN EVEN NUMBER

!     INTEGER, PARAMETER :: MAXSTREAMS_BRDF = 2
      INTEGER, PARAMETER :: MAXSTREAMS_BRDF = 100    ! best

!  Maximum numbers for the MSR quadratures

      INTEGER, PARAMETER :: max_msrs_muquad  = 50
      INTEGER, PARAMETER :: max_msrs_phiquad = 100

!  Number of quadrature streams for internal WSA/BSA scaling
!    New, Version 2.7. User does not need to know this value.

      INTEGER, parameter :: MAXSTREAMS_SCALING = 24

!  Weighting functions
!  -------------------

!  Maximum number of profile/column weighting functions

      INTEGER, PARAMETER :: MAX_ATMOSWFS = 3

!  Maximum number of surface property weighting functions

      INTEGER, PARAMETER :: MAX_SURFACEWFS = 6

!  Maximum number of surface-leaving weighting functions

      INTEGER, PARAMETER :: MAX_SLEAVEWFS = 1

!  Maximum number of error messages

      INTEGER, PARAMETER :: MAX_MESSAGES = 25

!  Derived dimensions
!  ==================

!  Copy Beam dimensioning

      INTEGER, PARAMETER :: MAXBEAMS = MAX_SZANGLES

!  Copy viewing zenith angle dimension

      INTEGER, PARAMETER :: MAX_USER_STREAMS = MAX_USER_VZANGLES

!  Maximum possible geometries

      INTEGER, PARAMETER :: MAX_GEOMETRIES = &
                            MAX_USER_VZANGLES*MAX_USER_RELAZMS*MAX_SZANGLES

!  All streams

      INTEGER, PARAMETER :: MAX_ALLSTRMS = MAX_USER_STREAMS + MAXSTREAMS

!  All streams for the Legendre PI-matrix setup.
!   Refractive Goemetry setting: Watch out for Kill. Memory Hog
!      INTEGER, PARAMETER :: MAX_ALLSTRMS_P1 = &
!                            MAX_ALLSTRMS + MAXBEAMS*MAXLAYERS

!  All streams for the Legendre PI-matrix setup.
!   Straightline setting: This setting should avoid dimensioning error

      INTEGER, PARAMETER :: MAX_ALLSTRMS_P1 = MAX_ALLSTRMS + MAXBEAMS

!  Maximum number of moments in the diffuse field calculation
!   This is always 2*MAXSTREAMS, in case we need DELTA-M

      INTEGER, PARAMETER :: MAXMOMENTS = 2*MAXSTREAMS

!  Maximum number of Fourier components = 2*MAXSTREAMS - 1

      INTEGER, PARAMETER :: MAXFOURIER = 2*MAXSTREAMS - 1

!  Number of Stokes streams squared

      INTEGER, PARAMETER :: MAXSTOKES_SQ = 16

!  Half the number of BRDF azimuth quadratures

      INTEGER, PARAMETER :: MAXSTHALF_BRDF = MAXSTREAMS_BRDF / 2

!  Other derived dimensions

      INTEGER, PARAMETER :: MAXSTREAMS_2   = 2*MAXSTREAMS
      INTEGER, PARAMETER :: MAXSTREAMS_21  = 2*MAXSTREAMS - 1
      INTEGER, PARAMETER :: MAXSTREAMS_P1  = MAXSTREAMS + 1
      INTEGER, PARAMETER :: MAXSTREAMS_P2  = MAXSTREAMS + 2

      INTEGER, PARAMETER :: MAXSTRMSTKS    = MAXSTREAMS * MAXSTOKES
      INTEGER, PARAMETER :: MAXSTRMSTKS_2  = 2*MAXSTRMSTKS
      INTEGER, PARAMETER :: MAXSTRMSTKS_21 = 2*MAXSTRMSTKS - 1

      INTEGER, PARAMETER :: MAXSTRMSTKS_P1 = MAXSTRMSTKS + 1
      INTEGER, PARAMETER :: MAXSTRMSTKS_P2 = MAXSTRMSTKS + 2
      INTEGER, PARAMETER :: MAXSTRMSTKS_P4 = MAXSTRMSTKS + 4

      INTEGER, PARAMETER :: MAX_USTRMSTKS = MAX_USER_STREAMS * MAXSTOKES

!  Maximum number of eigenvalues

      INTEGER, PARAMETER :: MAXEVALUES = MAXSTRMSTKS

!  For the BVP problem

      INTEGER, PARAMETER :: MAXTOTAL = MAXLAYERS*MAXSTRMSTKS_2
      INTEGER, PARAMETER :: MAXBANDTOTAL = 9*MAXSTRMSTKS - 2

!  Not so far used

      INTEGER, PARAMETER :: MAX_PSOLS = 2
      INTEGER, PARAMETER :: MAX_SCATPSOLS = MAX_PSOLS

!  Format constants
!  ================

      CHARACTER (LEN=*), PARAMETER :: &
        FMT_HEADING = '( / T6, ''-----> '', A, /)'

      CHARACTER (LEN=*), PARAMETER :: &
        FMT_INTEGER = '(T6, A, T58, I10)'

      CHARACTER (LEN=*), PARAMETER :: &
        FMT_REAL    = '(T6, A, T58, 1PG14.6)'

      CHARACTER (LEN=*), PARAMETER :: &
        FMT_CHAR    = '(T6, A, T48, A20)'

      CHARACTER (LEN=*), PARAMETER :: &
        FMT_SECTION = '( / T6, ''-----> '', A, /)'
!     $                  '( // T6
!     $                  ''----------------------------------------'',
!     $                  ''-----------------------------------'',
!     $                      / T6 A,
!     $                      / T6
!     $                  ''----------------------------------------'',
!     $                  ''-----------------------------------'',
!     $                          / )' )

!  Debug write format (DWF) constants
!  ==================================

      CHARACTER (LEN=*), PARAMETER :: DWFL  = '(A,L1)'
      CHARACTER (LEN=*), PARAMETER :: DWFL1 = '(A,I3,A,L1)'
      CHARACTER (LEN=*), PARAMETER :: DWFL2 = '(2(A,I3),A,L1)'

      CHARACTER (LEN=*), PARAMETER :: DWFI  = '(A,I5)'
      CHARACTER (LEN=*), PARAMETER :: DWFI1 = '(A,I3,A,I5)'
      CHARACTER (LEN=*), PARAMETER :: DWFI2 = '(2(A,I3),A,I5)'

      CHARACTER (LEN=*), PARAMETER :: DWFR  = '(A,ES13.6E2)'
      CHARACTER (LEN=*), PARAMETER :: DWFR1 = '(A,I3,A,ES13.6E2)'
      CHARACTER (LEN=*), PARAMETER :: DWFR2 = '(2(A,I3),A,ES13.6E2)'
      CHARACTER (LEN=*), PARAMETER :: DWFR3 = '(3(A,I3),A,ES13.6E2)'
      CHARACTER (LEN=*), PARAMETER :: DWFR4 = '(4(A,I3),A,ES13.6E2)'
      CHARACTER (LEN=*), PARAMETER :: DWFR5 = '(5(A,I3),A,ES13.6E2)'
      CHARACTER (LEN=*), PARAMETER :: DWFR6 = '(6(A,I3),A,ES13.6E2)'
      CHARACTER (LEN=*), PARAMETER :: DWFR7 = '(7(A,I3),A,ES13.6E2)'

      CHARACTER (LEN=*), PARAMETER :: DWFR1_3 = '(A,I3,3(A,ES13.6E2))'

      CHARACTER (LEN=*), PARAMETER :: DWFC  = '(2A)'
      CHARACTER (LEN=*), PARAMETER :: DWFC1 = '(A,I3,2A)'
      CHARACTER (LEN=*), PARAMETER :: DWFC2 = '(2(A,I3),2A)'

!  numbers
!  =======

      DOUBLE PRECISION, PARAMETER :: &
        ONE = 1.0D0, ZERO = 0.0D0,  ONEP5 = 1.5D0
      DOUBLE PRECISION, PARAMETER :: &
        TWO = 2.0D0, THREE = 3.0D0, FOUR = 4.0D0
      DOUBLE PRECISION, PARAMETER :: &
        QUARTER = 0.25D0, HALF = 0.5D0
      DOUBLE PRECISION, PARAMETER :: &
        MINUS_ONE = -ONE
      DOUBLE PRECISION, PARAMETER :: &
        MINUS_TWO = -TWO
      DOUBLE PRECISION, PARAMETER :: &
        DEG_TO_RAD = 1.7453292519943D-02
      DOUBLE PRECISION, PARAMETER :: &
        PIE = 180.0D0*DEG_TO_RAD
      DOUBLE PRECISION, PARAMETER :: &
        PI2 = 2.0D0 * PIE
      DOUBLE PRECISION, PARAMETER :: &
        PI4 = 4.0D0 * PIE
      DOUBLE PRECISION, PARAMETER :: &
        PIO2 = HALF * PIE
      DOUBLE PRECISION, PARAMETER :: &
        PIO4 = QUARTER * PIE
      DOUBLE PRECISION, PARAMETER :: &
        EPS3 = 0.001D0
      DOUBLE PRECISION, PARAMETER :: &
        EPS4 = 0.0001D0
      DOUBLE PRECISION, PARAMETER :: &
        EPS5 = 0.00001D0
      DOUBLE PRECISION, PARAMETER :: &
        SMALLNUM = 1.0D-9 !1.0D-15
      DOUBLE PRECISION, PARAMETER :: &
        BIGEXP = 32.0D0

!  Rob fix 5/6/13 - Taylor series limiting values

      DOUBLE PRECISION, PARAMETER :: TAYLOR_SMALL = 1.0D-04
      DOUBLE PRECISION, PARAMETER :: TAYLOR_LARGE = 1.0d+04

!  Control for Using L'Hopital's Rule
!   Changed, January 2009 for Version 2.4..........

      DOUBLE PRECISION, PARAMETER :: HOPITAL_TOLERANCE = EPS3
!      DOUBLE PRECISION, PARAMETER :: HOPITAL_TOLERANCE = EPS5

!  Control for limits of single scatter albedo

      DOUBLE PRECISION, PARAMETER :: OMEGA_SMALLNUM = 1.0D-15

!  Control for limits of extinction optical depth along solar path

      DOUBLE PRECISION, PARAMETER :: MAX_TAU_SPATH = 32.0D0

!  Control for limits of extinction optical depth along USER paths

      DOUBLE PRECISION, PARAMETER :: MAX_TAU_UPATH = 32.0D0

!  Control for limits of extinction optical depth along QUADRATURE paths

      DOUBLE PRECISION, PARAMETER :: MAX_TAU_QPATH = 32.0D0

!  error indices
!  =============

      INTEGER, PARAMETER :: VLIDORT_SERIOUS  = 4
      INTEGER, PARAMETER :: VLIDORT_WARNING  = 3
      INTEGER, PARAMETER :: VLIDORT_INFO     = 2
      INTEGER, PARAMETER :: VLIDORT_DEBUG    = 1
      INTEGER, PARAMETER :: VLIDORT_SUCCESS  = 0

!  directional indices

      INTEGER, PARAMETER :: UPIDX  = 1
      INTEGER, PARAMETER :: DNIDX  = 2

!  surface indices
!  ---------------

!  These refer to the BRDF kernel functions currently included.
!  Rob Extension 12/2/14. BPDF Kernels (replaces BPDF2009)

      INTEGER, PARAMETER :: LAMBERTIAN_IDX       = 1
      INTEGER, PARAMETER :: ROSSTHIN_IDX         = 2
      INTEGER, PARAMETER :: ROSSTHICK_IDX        = 3
      INTEGER, PARAMETER :: LISPARSE_IDX         = 4
      INTEGER, PARAMETER :: LIDENSE_IDX          = 5
      INTEGER, PARAMETER :: HAPKE_IDX            = 6
      INTEGER, PARAMETER :: ROUJEAN_IDX          = 7
      INTEGER, PARAMETER :: RAHMAN_IDX           = 8
      INTEGER, PARAMETER :: COXMUNK_IDX          = 9
      INTEGER, PARAMETER :: GISSCOXMUNK_IDX      = 10
      INTEGER, PARAMETER :: GISSCOXMUNK_CRI_IDX  = 11

!  New for Version 2.4RTC and up to Version 2.6.
!      INTEGER, PARAMETER :: BPDF2009_IDX         = 12

!  New BPDF functions for Version 2.7
!   These match the options in LIDORT Version 3.7.

      INTEGER, PARAMETER :: BPDFSOIL_IDX         = 12
      INTEGER, PARAMETER :: BPDFVEGN_IDX         = 13
      INTEGER, PARAMETER :: BPDFNDVI_IDX         = 14

!  New Cox-Munk functions for Version 2.7

      INTEGER, PARAMETER :: NewCMGLINT_IDX       = 15
      INTEGER, PARAMETER :: NewGCMGLINT_IDX      = 16

!  New for Version 2.8. Introduced 22 February 2016. R. Spurr
!  ----------------------------------------------------------

!  This is the Old Ross-Thick kernel with a hot-spot modification
!     Derives from the following references.

!  F. M. Breon, F. Maignan, M. Leroy and I. Grant, 
!    "Analysis of hot spot directional siganatures measured from space",
!      J. Geophys. Res., 107, D16, 4282, (2002)

!  E. Vermote C. Justice, and F. M. Breon,
!    "Towards a generalized approach for correction of the BRDF effect in MODIS reflectances"
!      IEEE Trans. Geo. Rem. Sens., 10.1109/TGRS.2008.2005997 (2008)
!   -- Ross-Thick Hotspot Kernel from

      INTEGER, PARAMETER :: RTKHOTSPOT_IDX       = 17

!  This is a "Modified Fresnel" kernel developed for polarized reflectances
!   Taken from the following reference.

!  P. Litvinov, O. Hasekamp and B. Cairns,
!    "Models for surface reflection of radiance and polarized radiance: Comparison
!     with airborne multi-angle photopolarimetric measurements and implications for
!     modeling top-of-atmopshere measurements,
!       Rem. Sens. Env., 115, 781-792, (2011).

      INTEGER, PARAMETER :: MODFRESNEL_IDX       = 18

!  1/31/21, Version 2.8.3. Analytical Model for Snow BRDF.
!     -- Kokhanovsky and Breon, IEEE GeoScience & Remote Sensing Letters, Vol 9(5), 928-932 (2012)
!     -- First introduced to VLIDORT, 18 November 2020.

      INTEGER, PARAMETER :: SNOWBRDF_IDX         = 19

!  1/31/21, Version 2.8.3. Maximum index
!    -- Revised to take into account New snow brdf index.

      INTEGER, PARAMETER :: MAXBRDF_IDX = SNOWBRDF_IDX

!  End of file.

      END MODULE vlidort_pars_m
