module TWOSTREAM_SurfaceMod

   USE TWOSTREAM_Mod

   implicit NONE

   integer, parameter :: LAMBERTIAN = 1
   integer, parameter :: GISSCOXMNK = 2
   integer, parameter :: MODIS_ROSSTHICK_LISPARSE = 3


   PUBLIC  TWOSTREAM_SurfaceLamb
   PUBLIC  TWOSTREAM_LANDMODIS

   TYPE TWOSTREAM_Surface

      type(TWOSTREAM)               :: Base
      integer                       :: sfc_type = -1
      real*8                        :: albedo

      real*8                        :: solar_zenith
      real*8                        :: relat_azimuth
      real*8                        :: sensor_zenith

   END TYPE TWOSTREAM_Surface

   contains

   Subroutine TWOSTREAM_SurfaceLamb(self, albedo, solar_zenith, &
                                  sensor_zenith, relative_azimuth)

         type(TWOSTREAM_Surface),intent(inout)  :: self
         real*8,  intent(in)  :: albedo
         real*8,  intent(in)  :: solar_zenith
         real*8,  intent(in)  :: sensor_zenith
         real*8,  intent(in)  :: relative_azimuth


         self%sfc_type      = 1
         self%albedo        = albedo
         self%solar_zenith  = solar_zenith
         self%sensor_zenith = sensor_zenith
         self%relat_azimuth = relative_azimuth

   end subroutine TWOSTREAM_SurfaceLamb

!.................................................................................

   Subroutine TWOSTREAM_LANDMODIS(self, solar_zenith, sensor_zenith, relative_azimuth, &
                                fiso,fgeo,fvol, param,rc)

      USE TWOSTREAM_PARS
      USE twostream_brdf_supplement_m

      implicit NONE
      type(LIDORT_Surface), intent(inout)   :: self
      real*8, intent(in)                    :: solar_zenith
      real*8, intent(in)                    :: sensor_zenith
      real*8, intent(in)                    :: relative_azimuth
      real*8, intent(in)                    :: fiso
      real*8, intent(in)                    :: fgeo
      real*8, intent(in)                    :: fvol
      real*8, intent(in), dimension(:)      :: param
      integer, intent(out)                  :: rc     ! error code

!     Local variables
!     ---------------
      real*8, dimension(MAX_USER_OBSGEOMS,3)             :: USER_OBSGEOMS
      real*8, dimension(MAXBEAMS)                        :: BEAM_SZAS
      real*8, dimension(MAX_USER_ANGLES)                 :: USER_ANGLES
      real*8, dimension(MAX_USER_RELAZMS)                :: USER_RELAZMS

!      BRDF variables
!     ---------------

      logical                               :: DO_SHADOW_EFFECT
      logical,dimension(MAX_BRDF_KERNELS)   ::   LAMBERTIAN_KERNEL_FLAG

      integer                                       ::   N_BRDF_KERNELS
      integer                                       ::   NSTREAMS_BRDF
      character(len=10),dimension(MAX_BRDF_KERNELS) ::   BRDF_NAMES
      integer,dimension(MAX_BRDF_KERNELS)           ::   WHICH_BRDF
      double precision, dimension(MAX_BRDF_KERNELS) ::   BRDF_FACTORS
      integer, dimension(MAX_BRDF_KERNELS)          ::   N_BRDF_PARAMETERS
      double precision, dimension(MAX_BRDF_KERNELS, MAX_BRDF_PARAMETERS)  ::   BRDF_PARAMETERS

      INTEGER       :: STATUS_BRDFSUP
      CHARACTER*100 :: MESSAGE_BRDF, ACTION_BRDF

      rc = 0

      self%sfc_type      = 3
      self%solar_zenith  = solar_zenith
      self%sensor_zenith = sensor_zenith
      self%relat_azimuth = relative_azimuth

      USER_OBSGEOMS = 0.0
      USER_OBSGEOMS(1,1) = solar_zenith
      USER_OBSGEOMS(1,2) = sensor_zenith
      USER_OBSGEOMS(1,3) = relat_azimuth
      BEAM_SZAS = 0.0
      BEAM_SZAS(1) = solar_zenith
      USER_ANGLES = 0.0
      USER_ANGLES(1) = sensor_zenith
      USER_RELAZMS = 0.0
      USER_RELAZMS(1) = relat_azimuth



      self%Base%TSIO%DO_BRDF_SURFACE = .true.
      DO_SHADOW_EFFECT = .false.               ! For Cox-Munk
      N_BRDF_KERNELS = 3                       ! Number of BRDF kernels (max 3)
      NSTREAMS_BRDF = 50                       ! Number of azimuth quadrature streams for BRDF
    
      ! For each BRDF_KERNELS specify the name, factor and parameter
      !-------------------------------------------------------------
      BRDF_NAMES(1) = 'Lambertian'   ! 0 free parameter
      WHICH_BRDF(1) =  1             ! specified in LIDORT.pars.f90
      BRDF_FACTORS(1) =  fiso        ! From MODIS MOD43
      N_BRDF_PARAMETERS(1) = 0
      BRDF_PARAMETERS(1,1) = 0.0
      BRDF_PARAMETERS(1,2) = 0.0
      BRDF_PARAMETERS(1,3) = 0.0
      LAMBERTIAN_KERNEL_FLAG(1) = .true. ! set .true. if BRDF_NAME(I) = 'Lambertian'

      BRDF_NAMES(2) = 'Ross-thick'   ! 0 free parameter
      WHICH_BRDF(2) =  3             ! specified in LIDORT.pars.f90
      BRDF_FACTORS(2) =  fvol        ! From MODIS
      N_BRDF_PARAMETERS(2) = 0
      BRDF_PARAMETERS(2,1) = 0.0
      BRDF_PARAMETERS(2,2) = 0.0
      BRDF_PARAMETERS(2,3) = 0.0
      LAMBERTIAN_KERNEL_FLAG(2) = .false.

      BRDF_NAMES(3) = 'Li-sparse'   ! 2 free parameters
      WHICH_BRDF(3) =  4            ! specified in LIDORT.pars.f90
      BRDF_FACTORS(3) =  fgeo       ! From MODIS
      N_BRDF_PARAMETERS(3) = 2
      BRDF_PARAMETERS(3,1) = param(1)    ! h/b (relative height) -> (h is the height-to-center
                                         ! of the spheroids from the ground, b is the vertical
                                         !  and r the horizontal spheroid radius
      BRDF_PARAMETERS(3,2) = param(2)    ! b/r (crown shape) (MODIS values and Wanner et al., 1995 Sect 5.1)
      BRDF_PARAMETERS(3,3) = 0.0
      LAMBERTIAN_KERNEL_FLAG(3) = .false.


      !Call BRDF supplement
      CALL TWOSTREAM_BRDFMASTER &
            ( MAXBEAMS, MAX_USER_ANGLES, MAX_USER_OBSGEOMS, & ! Dimensions !@@
              MAXSTREAMS_BRDF, MAX_BRDF_KERNELS,            & ! Dimensions
              MAX_BRDF_PARAMETERS,                          & ! Dimensions
              self%Base%TSIO%DO_SOLAR_SOURCES,              &
              self%Base%TSIO%DO_USER_OBSGEOMS,              & ! Inputs !@@
              LAMBERTIAN_KERNEL_FLAG,                       & ! Inputs
              DO_SHADOW_EFFECT,                             & 
              self%Base%VSIO%DO_SURFACE_EMISSION,           & ! Inputs
              self%Base%NBEAMS, self%Base%N_USER_ANGLES,    &
              self%Base%N_USER_OBSGEOMS,                    & ! Inputs !@@
              BEAM_SZAS, USER_ANGLES, USER_OBSGEOMS,        & ! Inputs !@@
              self%Base%TSIO%STREAM_VALUE, NSTREAMS_BRDF,   & ! Inputs
              N_BRDF_KERNELS, WHICH_BRDF, BRDF_FACTORS,     & ! Inputs
              N_BRDF_PARAMETERS, BRDF_PARAMETERS,           & ! Inputs
              self%Base%TSIO%BRDF_F_0,                      & ! Outputs
              self%Base%TSIO%BRDF_F,                        & ! Outputs
              self%Base%TSIO%UBRDF_F,                       & ! Outputs
              self%Base%TSIO%EMISSIVITY,                    & ! Outputs
              STATUS_BRDFSUP, MESSAGE_BRDF, ACTION_BRDF )     ! Outputs

      if ( STATUS_BRDFSUP .ne. 0 ) then
         write(*,*)'BRDF supplement Check failed'
         write(*,*)' - Print 1 Message and 1 Action'
         write(*,'(A)') TRIM(MESSAGE_BRDF)
         write(*,'(A)') TRIM(ACTION_BRDF)
         rc = 1
      endif

   end subroutine TWOSTREAM_LANDMODIS

   end module TWOSTREAM_SurfaceMod
