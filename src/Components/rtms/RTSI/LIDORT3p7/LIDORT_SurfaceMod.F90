module LIDORT_SurfaceMod

   USE LIDORT_Mod

   implicit NONE

   integer, parameter :: LAMBERTIAN = 1
   integer, parameter :: GISSCOXMNK = 2
   integer, parameter :: MODIS_ROSSTHICK_LISPARSE = 3


   PUBLIC  LIDORT_SurfaceLamb
   PUBLIC  LIDORT_GissCoxMunk
   PUBLIC  LIDORT_LANDMODIS
       
   TYPE LIDORT_Surface

      type(LIDORT)                  :: Base 
      integer                       :: sfc_type = -1
      real*8                        :: albedo 

      logical                       :: scalar = .false.
      
      real*8                        :: solar_zenith
      real*8                        :: relat_azimuth 
      real*8                        :: sensor_zenith 


   END TYPE LIDORT_Surface

   contains

   Subroutine LIDORT_SurfaceLamb(self, albedo, solar_zenith, &
                                  sensor_zenith, relative_azimuth)
         
         type(LIDORT_Surface),intent(inout)  :: self
         real*8,  intent(in)  :: albedo
         real*8,  intent(in)  :: solar_zenith
         real*8,  intent(in)  :: sensor_zenith
         real*8,  intent(in)  :: relative_azimuth
         
       
         self%sfc_type      = 1
         self%albedo        = albedo
         self%solar_zenith  = solar_zenith
         self%sensor_zenith = sensor_zenith
         self%relat_azimuth = relative_azimuth
      
   end subroutine LIDORT_SurfaceLamb
!.........................................................................
   Subroutine LIDORT_GissCoxMunk(self,U10m,V10m,m, solar_zenith, &
                                  sensor_zenith, relative_azimuth, BRDF,rc) 
 
      USE LIDORT_PARS
      USE LIDORT_AUX
      USE BRDF_SUP_MOD

      implicit NONE

      type(LIDORT_Surface), intent(inout)    :: self
      real*8, intent(in)                     :: solar_zenith
      real*8, intent(in)                     :: sensor_zenith
      real*8, intent(in)                     :: relative_azimuth 
      real*8, intent(in)                     :: U10m
      real*8, intent(in)                     :: V10m
      real*8, intent(in)                     :: m
      integer, intent(out)                   :: rc     ! error code
      real*8, intent(out)                    :: BRDF     
    
      !Input parameters BRDF supplement LIDORT:
      TYPE LIDORT_BRDF
       !  BRDF input structures
         TYPE(BRDF_Sup_Inputs)                :: BRDF_Sup_In
         TYPE(BRDF_Input_Exception_Handling)  :: BRDF_Sup_InputStatus
      
       !  LIDORT BRDF output of the supp but input for LIDORT scat
       !  TYPE(BRDF_Sup_Outputs)               :: BRDF_Sup_Out
      END TYPE LIDORT_BRDF

      type(LIDORT_BRDF)   :: LBRDF

!     Local variables 
!     ---------------
!     BRDF variables
!     ---------------
      logical                                         ::   DO_BRDF_SURFACE
      logical                                         ::   DO_USER_STREAMS
      logical,dimension(MAX_BRDF_KERNELS)             ::   LAMBERTIAN_KERNEL_FLAG
      logical                                         ::   DO_SURFACE_EMISSION 
      logical                                         ::   DO_SHADOW_EFFECT    ! Only for Cox-Munk type kernels
      logical                                         ::   DO_COXMUNK_DBMS     ! Only for Cox-Munk type kernels
      integer                                         ::   N_BRDF_KERNELS 
      integer                                         ::   NSTREAMS_BRDF
      character(len=10),dimension(MAX_BRDF_KERNELS)   ::   BRDF_NAMES
      integer,dimension(MAX_BRDF_KERNELS)             ::   WHICH_BRDF
      double precision, dimension(MAX_BRDF_KERNELS)   ::   BRDF_FACTORS
      integer, dimension(MAX_BRDF_KERNELS)            ::   N_BRDF_PARAMETERS
      double precision, dimension(MAX_BRDF_KERNELS, MAX_BRDF_PARAMETERS)    ::   BRDF_PARAMETERS
      logical                                         ::   DO_DEBUG_RESTORATION
      logical                                         ::   DO_EXACTONLY
      integer                                         ::   N_MOMENTS_INPUT
      logical                                         ::   DO_MSRCORR
      logical                                         ::   DO_MSRCORR_DBONLY  
      integer                                         ::   GLITTER_MSRCORR_ORDER
      integer                                         ::   GLITTER_MSRCORR_NMUQUAD
      integer                                         ::   GLITTER_MSRCORR_NPHIQUAD

      double precision                                ::  sigma2
!      double precision                               ::  m
      double precision                                ::  W10m

      rc = 0

      self%sfc_type      = 2
      self%solar_zenith  = solar_zenith
      self%sensor_zenith = sensor_zenith
      self%relat_azimuth = relative_azimuth               

      DO_USER_STREAMS     = .true.       ! Use user-defined viewing zenith angles
      DO_BRDF_SURFACE     = .true. 
      DO_SURFACE_EMISSION = .false.  ! Calculation of surface thermal emission

      N_BRDF_KERNELS            = 1             ! Number of BRDF kernels (max 3)
      NSTREAMS_BRDF             = 100             ! Number of BRDF azimuth angles
      LAMBERTIAN_KERNEL_FLAG(1) = .false. ! set .true. if BRDF_NAME(I) = 'Lambertian'

      ! Calculation of the slope-squared variance (fct of W : wind speed (m/s))
      ! ---------------------------------------------------------
      W10m = sqrt(U10m*U10m+V10m*V10m)
      if (W10m<0)   W10m= 0
      if (W10m>50)  W10m= 50
      sigma2 = 0.003 + 0.00512 * W10m
!      print*, 'sigma2', sigma2

      ! air-water refractive index m
      ! -------------------
!      m = 1.33   

      ! For each BRDF_KERNELS specify the name, factor and parameter  
      !-------------------------------------------------------------   
      BRDF_NAMES(1) = 'GissCoxMnk' 
!      BRDF_NAMES(1) ='Cox-Munk  '  
      WHICH_BRDF(1) =  10            ! specified in LIDORT.pars.f90
!      WHICH_BRDF(1) =  9   
      BRDF_FACTORS(1)      =  1.   
      N_BRDF_PARAMETERS(1) = 3
      BRDF_PARAMETERS(1,1) = 0.5*sigma2
      BRDF_PARAMETERS(1,2) = m
      BRDF_PARAMETERS(1,3) = 0

      ! Only for Cox-Munk type kernels
      DO_SHADOW_EFFECT = .true.      ! Shadow effect for glitter kernels
!      DO_DBONLY = .false.        ! only direct bounce BRDF (no Fourier terms calc)

      DO_MSRCORR         = .true.  ! Do multiple reflectance correction for glitter kernels
      DO_MSRCORR_DBONLY  = .false. ! Do multiple reflectance correction 
                                               !for only the exact term of glitter kernels

      GLITTER_MSRCORR_ORDER    = 0   !Order of correction for multiple reflectance computations
                                    !(0 = no corr), 1,,2,3 ...
      GLITTER_MSRCORR_NMUQUAD  = 40    ! Number of angles used in zenith integration
                                    ! during  multiple corr
      GLITTER_MSRCORR_NPHIQUAD = 100   ! Number of angles used in azimuthal integration
                                    ! during  multiple corr

      DO_USER_STREAMS          = .true.       ! Use user-defined viewing zenith angles

                                           
      ! Copy in data structure
      ! ---------------------

      LBRDF%BRDF_Sup_In%BS_DO_USER_STREAMS     = DO_USER_STREAMS
      LBRDF%BRDF_Sup_In%BS_DO_BRDF_SURFACE     = DO_BRDF_SURFACE
      LBRDF%BRDF_Sup_In%BS_DO_SURFACE_EMISSION = DO_SURFACE_EMISSION

      ! Angles

      LBRDF%BRDF_Sup_In%BS_NSTREAMS             = self%Base%NSTREAMS
      LBRDF%BRDF_Sup_In%BS_NBEAMS               = self%Base%NBEAMS
      LBRDF%BRDF_Sup_In%BS_BEAM_SZAS(1)         = solar_zenith
      LBRDF%BRDF_Sup_In%BS_N_USER_RELAZMS       = self%Base%N_USER_RELAZMS
      LBRDF%BRDF_Sup_In%BS_USER_RELAZMS(1)      = relative_azimuth
      LBRDF%BRDF_Sup_In%BS_N_USER_STREAMS       = self%Base%N_USER_STREAMS
      LBRDF%BRDF_Sup_In%BS_USER_ANGLES_INPUT(1) = sensor_zenith

      ! BRDF inputs

      LBRDF%BRDF_Sup_In%BS_N_BRDF_KERNELS         = N_BRDF_KERNELS
      LBRDF%BRDF_Sup_In%BS_BRDF_NAMES             = BRDF_NAMES
      LBRDF%BRDF_Sup_In%BS_WHICH_BRDF             = WHICH_BRDF
      LBRDF%BRDF_Sup_In%BS_N_BRDF_PARAMETERS      = N_BRDF_PARAMETERS
      LBRDF%BRDF_Sup_In%BS_BRDF_PARAMETERS        = BRDF_PARAMETERS
      LBRDF%BRDF_Sup_In%BS_LAMBERTIAN_KERNEL_FLAG = LAMBERTIAN_KERNEL_FLAG
      LBRDF%BRDF_Sup_In%BS_BRDF_FACTORS           = BRDF_FACTORS
      LBRDF%BRDF_Sup_In%BS_NSTREAMS_BRDF          = NSTREAMS_BRDF
      LBRDF%BRDF_Sup_In%BS_DO_SHADOW_EFFECT       = DO_SHADOW_EFFECT
!      LBRDF%BRDF_Sup_In%BS_DO_EXACTONLY           = DO_EXACTONLY
     

      LBRDF%BRDF_Sup_In%BS_DO_GLITTER_MSRCORR           = DO_MSRCORR
      LBRDF%BRDF_Sup_In%BS_DO_GLITTER_MSRCORR_DBONLY = DO_MSRCORR_DBONLY
      LBRDF%BRDF_Sup_In%BS_GLITTER_MSRCORR_ORDER        = GLITTER_MSRCORR_ORDER
      LBRDF%BRDF_Sup_In%BS_GLITTER_MSRCORR_NMUQUAD      = GLITTER_MSRCORR_NMUQUAD
      LBRDF%BRDF_Sup_In%BS_GLITTER_MSRCORR_NPHIQUAD     = GLITTER_MSRCORR_NPHIQUAD


      if ( N_BRDF_KERNELS .GT. MAX_BRDF_KERNELS )       rc = 1
      if ( NSTREAMS_BRDF .GT. MAXSTREAMS_BRDF )         rc = 2
      if ( WHICH_BRDF(1) .GT. MAXBRDF_IDX )             rc = 4


     
      DO_DEBUG_RESTORATION = .false.
      N_MOMENTS_INPUT = 2 * LBRDF%BRDF_Sup_In%BS_NSTREAMS - 1

      call BRDF_MAINMASTER (DO_DEBUG_RESTORATION, &         ! Inputs
                             N_MOMENTS_INPUT, &              ! Inputs
                             LBRDF%BRDF_Sup_In, &           ! Inputs
                             self%Base%VIO%BRDF_Sup_Out, &  ! Outputs
                           self%Base%VIO%BRDF_Sup_OutputStatus)          ! Outputs
       BRDF = self%Base%VIO%BRDF_Sup_Out%BS_DBOUNCE_BRDFUNC(1,1,1)     
       print*, 'ok BRDF', self%Base%VIO%BRDF_Sup_Out%BS_DBOUNCE_BRDFUNC(1,1,1)       
   end subroutine LIDORT_GissCoxMunk

!.................................................................................

   Subroutine LIDORT_LANDMODIS(self, solar_zenith, sensor_zenith, relative_azimuth, &
                                fiso,fgeo,fvol, param,rc) 
 
      USE LIDORT_PARS
      USE LIDORT_AUX
      USE BRDF_SUP_MOD

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
!      real*8, intent(out)                   :: BRDF     

      !Input parameters for LIDORT BRDF supplement:
      TYPE LIDORT_BRDF
       !  BRDF input structures
         TYPE(BRDF_Sup_Inputs)                :: BRDF_Sup_In
         TYPE(BRDF_Input_Exception_Handling)  :: BRDF_Sup_InputStatus
      
       !  LIDORT BRDF output of the supp but input for LIDORT scat
       !  TYPE(BRDF_Sup_Outputs)               :: BRDF_Sup_Out
      END TYPE LIDORT_BRDF

      type(LIDORT_BRDF)   :: LBRDF

       
!     Local variables 
!     ---------------
!      BRDF variables
!     ---------------
      logical                             ::   DO_BRDF_SURFACE
      logical                             ::   DO_USER_STREAMS
      logical,dimension(MAX_BRDF_KERNELS) ::   LAMBERTIAN_KERNEL_FLAG
      logical                             ::   DO_SURFACE_EMISSION 
      logical                             ::   DO_SOLAR_SOURCES
      logical                             ::   DO_USER_OBSGEOMS

      integer                                       ::   N_BRDF_KERNELS 
      integer                                       ::   NSTREAMS_BRDF
      character(len=10),dimension(MAX_BRDF_KERNELS) ::   BRDF_NAMES
      integer,dimension(MAX_BRDF_KERNELS)           ::   WHICH_BRDF
      double precision, dimension(MAX_BRDF_KERNELS) ::   BRDF_FACTORS
      integer, dimension(MAX_BRDF_KERNELS)          ::   N_BRDF_PARAMETERS
      double precision, dimension(MAX_BRDF_KERNELS, MAX_BRDF_PARAMETERS)  ::   BRDF_PARAMETERS
      logical                                       ::   DO_DEBUG_RESTORATION
      integer                                       ::   N_MOMENTS_INPUT
      integer                                       ::   i
      CHARACTER (LEN=120)                           :: MESSAGES ( 0:MAX_MESSAGES )
      CHARACTER (LEN=120)                           :: ACTIONS ( 0:MAX_MESSAGES )

      rc = 0

      self%sfc_type      = 3
      self%solar_zenith  = solar_zenith
      self%sensor_zenith = sensor_zenith
      self%relat_azimuth = relative_azimuth  
             

      DO_BRDF_SURFACE  = .true. 
      DO_USER_STREAMS  = .true.      ! Use user-defined viewing zenith angles
      DO_SOLAR_SOURCES = .true.      ! TRUE for sunlight, may be TRUE or FALSE in thermal regime
      DO_USER_OBSGEOMS = .true.      ! BRDF is only calculated for specific azimuth angle

      DO_SURFACE_EMISSION = .false.  ! Calculation of surface thermal emission


      N_BRDF_KERNELS = 3             ! Number of BRDF kernels (max 3)
      NSTREAMS_BRDF = 50             ! Number of azimuth quadrature streams for BRDF
                                     ! Not really needed for radiance calc if DO_USER_OBSGEOMS is TRUE 
          

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

      ! Copy in data structure
      ! ---------------------

      LBRDF%BRDF_Sup_In%BS_DO_USER_STREAMS     = DO_USER_STREAMS
      LBRDF%BRDF_Sup_In%BS_DO_BRDF_SURFACE     = DO_BRDF_SURFACE
      LBRDF%BRDF_Sup_In%BS_DO_SURFACE_EMISSION = DO_SURFACE_EMISSION
      LBRDF%BRDF_Sup_In%BS_DO_SOLAR_SOURCES    = DO_SOLAR_SOURCES   
      LBRDF%BRDF_Sup_In%BS_DO_USER_OBSGEOMS    = DO_USER_OBSGEOMS   

      ! Angles

      LBRDF%BRDF_Sup_In%BS_NSTREAMS             = self%Base%NSTREAMS
      LBRDF%BRDF_Sup_In%BS_NBEAMS               = self%Base%NBEAMS
      LBRDF%BRDF_Sup_In%BS_BEAM_SZAS(1)         = solar_zenith
      LBRDF%BRDF_Sup_In%BS_N_USER_RELAZMS       = self%Base%N_USER_RELAZMS
      LBRDF%BRDF_Sup_In%BS_USER_RELAZMS(1)      = relative_azimuth
      LBRDF%BRDF_Sup_In%BS_N_USER_STREAMS       = self%Base%N_USER_STREAMS
      LBRDF%BRDF_Sup_In%BS_USER_ANGLES_INPUT(1) = sensor_zenith
      LBRDF%BRDF_Sup_In%BS_N_USER_OBSGEOMS      = self%Base%N_USER_OBSGEOMS 
      LBRDF%BRDF_Sup_In%BS_USER_OBSGEOMS(1,1)   = solar_zenith
      LBRDF%BRDF_Sup_In%BS_USER_OBSGEOMS(1,2)   = sensor_zenith
      LBRDF%BRDF_Sup_In%BS_USER_OBSGEOMS(1,3)   = relative_azimuth                       
    
      ! BRDF inputs

      LBRDF%BRDF_Sup_In%BS_N_BRDF_KERNELS         = N_BRDF_KERNELS
      LBRDF%BRDF_Sup_In%BS_BRDF_NAMES             = BRDF_NAMES
      LBRDF%BRDF_Sup_In%BS_WHICH_BRDF             = WHICH_BRDF
      LBRDF%BRDF_Sup_In%BS_N_BRDF_PARAMETERS      = N_BRDF_PARAMETERS
      LBRDF%BRDF_Sup_In%BS_BRDF_PARAMETERS        = BRDF_PARAMETERS
      LBRDF%BRDF_Sup_In%BS_LAMBERTIAN_KERNEL_FLAG = LAMBERTIAN_KERNEL_FLAG
      LBRDF%BRDF_Sup_In%BS_BRDF_FACTORS           = BRDF_FACTORS
      LBRDF%BRDF_Sup_In%BS_NSTREAMS_BRDF          = NSTREAMS_BRDF


      LBRDF%BRDF_Sup_In%BS_DO_DIRECTBOUNCE_ONLY   = .false.

      ! The following is needed for Cox-munck type only
      ! set to initialization values
      !-----------------------------------------------
      LBRDF%BRDF_Sup_In%BS_DO_SHADOW_EFFECT          = .false.  

      LBRDF%BRDF_Sup_In%BS_DO_GLITTER_MSRCORR        = .false.
      LBRDF%BRDF_Sup_In%BS_DO_GLITTER_MSRCORR_DBONLY = .false. 
      LBRDF%BRDF_Sup_In%BS_GLITTER_MSRCORR_ORDER     = 0
      LBRDF%BRDF_Sup_In%BS_GLITTER_MSRCORR_NMUQUAD   = 0
      LBRDF%BRDF_Sup_In%BS_GLITTER_MSRCORR_NPHIQUAD  = 0

      LBRDF%BRDF_Sup_In%BS_DO_NewCMGLINT             = .false.
      LBRDF%BRDF_Sup_In%BS_SALINITY                  = 0
      LBRDF%BRDF_Sup_In%BS_WAVELENGTH                = 0

      LBRDF%BRDF_Sup_In%BS_WINDSPEED                 = 0
      LBRDF%BRDF_Sup_In%BS_WINDDIR                   = 0

      LBRDF%BRDF_Sup_In%BS_DO_GlintShadow            = .false.
      LBRDF%BRDF_Sup_In%BS_DO_FoamOption             = .false.
      LBRDF%BRDF_Sup_In%BS_DO_FacetIsotropy          = .false.

      ! WSA and BSA scaling options.
      ! WSA = White-sky albedo. BSA = Black-sky albedo.
      ! Not implemented for now.  Could be tested
      !--------------------------------------------
      LBRDF%BRDF_Sup_In%BS_DO_WSABSA_OUTPUT = .true.
      LBRDF%BRDF_Sup_In%BS_DO_WSA_SCALING   = .false.
      LBRDF%BRDF_Sup_In%BS_DO_BSA_SCALING   = .false.
      LBRDF%BRDF_Sup_In%BS_WSA_VALUE        = 0
      LBRDF%BRDF_Sup_In%BS_BSA_VALUE        = 0

      !  Exception handling
      ! ----------------------------
      MESSAGES(0)     = 'Successful Read of LIDORT Input file'
      ACTIONS(0)      = 'No Action required for this Task'

      LBRDF%BRDF_Sup_InputStatus%BS_STATUS_INPUTREAD = LIDORT_SUCCESS

      LBRDF%BRDF_Sup_InputStatus%BS_NINPUTMESSAGES   = 0
      LBRDF%BRDF_Sup_InputStatus%BS_INPUTMESSAGES    = MESSAGES
      LBRDF%BRDF_Sup_InputStatus%BS_INPUTACTIONS     = ACTIONS
      
     
      ! Do some checks to make sure parameters are not out of range
      !------------------------------------------------------------
      if ( N_BRDF_KERNELS .GT. MAX_BRDF_KERNELS )       rc = 1
      if ( NSTREAMS_BRDF .GT. MAXSTREAMS_BRDF )         rc = 2

      do i =1, N_BRDF_KERNELS
         if ( WHICH_BRDF(i) .GT. MAXBRDF_IDX )          rc = 4
      end do
 
      ! Debug flag for restoration
      DO_DEBUG_RESTORATION = .false.
      ! Number of moments (only used for restoration debug)
      N_MOMENTS_INPUT = 2 * LBRDF%BRDF_Sup_In%BS_NSTREAMS - 1

      ! Run the LIDORT Surface Module
      ! -------------------------------------------------
      call BRDF_MAINMASTER (DO_DEBUG_RESTORATION, &                ! Inputs
                             N_MOMENTS_INPUT, &                     ! Inputs
                             LBRDF%BRDF_Sup_In, &                  ! Inputs
                             self%Base%VIO%BRDF_Sup_Out, &         ! Outputs
                             self%Base%VIO%BRDF_Sup_OutputStatus)          ! Outputs
!       BRDF = self%Base%VIO%BRDF_Sup_Out%BS_DBOUNCE_BRDFUNC(1,1,1)     
!       print*, 'ok BRDF', self%Base%VIO%BRDF_Sup_Out%BS_DBOUNCE_BRDFUNC(1,1,1)             
   end subroutine LIDORT_LANDMODIS

   end module LIDORT_SurfaceMod
