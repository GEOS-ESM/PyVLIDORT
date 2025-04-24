      module TWOSTREAM_ScatMod

      USE TWOSTREAM_PARS
      USE twostream_brdf_supplement_m
      USE twostream_master_m

      USE TWOSTREAM_Mod
      USE TWOSTREAM_SurfaceMod

      implicit NONE

      PUBLIC  TWOSTREAM_Run

      type TWOSTREAM_scat

         real*8          :: wavelength          ! in [nm]
         integer         :: nMom                ! number of momemts read (phase function)
         real*8          :: MISSING             ! MISSING VALUE
         real*8, pointer :: rot(:)              ! rayleigh optical thickness
         real*8, pointer :: depol_ratio
         real*8, pointer :: alpha(:)            ! trace gas absorption optical thickness
         real*8, pointer :: tau(:)              ! aerosol tau
         real*8, pointer :: ssa(:)              ! aerosol ssa
         real*8, pointer ::   g(:)              ! aerosol asymmetry factor
         real*8, pointer ::  pe(:)              ! pressure    at layer edges [Pa]
         real*8, pointer ::  ze(:)              ! height      at layer edges [m]
         real*8, pointer ::  te(:)              ! temperature at layer edges [K]

         type(TWOSTREAM_Surface) :: Surface

      end type TWOSTREAM_scat

      type TWOSTREAM_output
         real*8     :: RADIANCE    ! TOA radiance
         real*8     :: REFLECTANCE ! TOA reflectance
      end type TWOSTREAM_output

      contains
!.............................................................................

      subroutine TWOSTREAM_Run (self, output, rc)
!
!     Computes radiances for a single wavelength, pixel. Optical properties
!     and met fields in self are assumed to have been updated with the
!     apropriate values.
!
      USE TWOSTREAM_PARS

      USE TWOSTREAM_MOD
      USE twostream_brdf_supplement_m

      USE twostream_master_m

      type(TWOSTREAM_scat),          intent(inout)   :: self        ! Contains most input
      type(TWOSTREAM_output),        intent(out)     :: output      ! contains output
      integer,                       intent(out)     :: rc


!     local variables
!     ---------------
      integer                                            :: i, j, k, l, m, n
      integer                                            :: NLAYERS
      real*8                                             :: ray_l
      real*8                                             :: tau_l
      real*8                                             :: ssa_l
      real*8                                             :: g_l
      real*8                                             :: tau_ext
      real*8                                             :: tau_scat
      real*8                                             :: ssa_tot
      real*8                                             :: raysmom2
      real*8                                             :: aerswt
      real*8                                             :: rayswt
      real*8                                             :: m1, m2

      real*8                                             :: LAMBERTIAN_ALBEDO

      real*8, dimension(0:MAXLAYERS)                     :: HEIGHT_GRID ( 0:MAXLAYERS )
      real*8, dimension(MAXLAYERS)                       :: assym_vert_input
      real*8, dimension(MAXLAYERS)                       :: deltau_vert_input
      real*8, dimension(MAXLAYERS)                       :: omega_total_input
      real*8, dimension(MAXLAYERS)                       :: d2s_scaling

      real*8, dimension(MAX_USER_OBSGEOMS,3)             :: USER_OBSGEOMS
      real*8, dimension(MAXBEAMS)                        :: BEAM_SZAS
      real*8, dimension(MAX_USER_ANGLES)                 :: USER_ANGLES
      real*8, dimension(MAX_USER_RELAZMS)                :: USER_RELAZMS

      real*8, dimension(MAX_GEOMETRIES)                  :: INTENSITY_TOA
      real*8, dimension(MAX_GEOMETRIES)                  :: INTENSITY_BOA 
      real*8, dimension(MAXBEAMS,2)                      :: FLUXES_TOA
      real*8, dimension(MAXBEAMS,2)                      :: FLUXES_BOA
      real*8, dimension(MAX_GEOMETRIES,0:MAXLAYERS)      :: RADLEVEL_UP 
      real*8, dimension(MAX_GEOMETRIES,0:MAXLAYERS)      :: RADLEVEL_DN 

!  Numbers (geometry)
!   N_GEOMETRIES = NBEAMS * N_USER_STREAMS * N_USER_RELAZMS (Lattice value)
!   N_GEOMETRIES = N_USER_OBSGEOMS                          (OBsGeom value)

      INTEGER       :: n_geometries


      real*8, parameter                                  :: pi = 4.*atan(1.0)

      logical                                            :: DO_UPWELLING          
      logical                                            :: DO_DNWELLING          
      logical                                            :: DO_PLANE_PARALLEL     
      logical                                            :: DO_2S_LEVELOUT         
      logical                                            :: DO_MVOUT_ONLY          
      logical                                            :: DO_ADDITIONAL_MVOUT    
      logical                                            :: DO_SOLAR_SOURCES       
      logical                                            :: DO_THERMAL_EMISSION    
      logical                                            :: DO_SURFACE_EMISSION    
      logical                                            :: DO_D2S_SCALING         
      logical                                            :: DO_USER_OBSGEOMS       
      logical                                            :: DO_SURFACE_LEAVING     
      logical                                            :: DO_SL_ISOTROPIC        
      logical                                            :: DO_PENTADIAG_INVERSE  
      logical                                            :: DO_BRDF_SURFACE 


      rc = 0

      if ( .not. self%Surface%Base%initialized ) then
        rc = 1
        return
      end if

      NLAYERS = self%Surface%Base%TSIO%NLAYERS

!                          Lambertian OR BRDF surface
!                          --------------------------

      if ( self%Surface%sfc_type == 1 ) then               ! Lambertian surface ?
         DO_BRDF_SURFACE = .false.
         LAMBERTIAN_ALBEDO = self%Surface%albedo           ! Lambertian (isotropic) input albedo
      else
         DO_BRDF_SURFACE = .true.
         LAMBERTIAN_ALBEDO = 0.0                           ! Use BRDF -> albedo set to 0
      end if


!                         Copy over Boolean Inputs
!                         ------------------------
      DO_UPWELLING           = self%Surface%Base%TSIO%DO_UPWELLING
      DO_DNWELLING           = self%Surface%Base%TSIO%DO_DNWELLING
      DO_PLANE_PARALLEL      = self%Surface%Base%TSIO%DO_PLANE_PARALLEL
      DO_2S_LEVELOUT         = self%Surface%Base%TSIO%DO_2S_LEVELOUT
      DO_MVOUT_ONLY          = self%Surface%Base%TSIO%DO_MVOUT_ONLY
      DO_ADDITIONAL_MVOUT    = self%Surface%Base%TSIO%DO_ADDITIONAL_MVOUT
      DO_SOLAR_SOURCES       = self%Surface%Base%TSIO%DO_SOLAR_SOURCES
      DO_THERMAL_EMISSION    = self%Surface%Base%TSIO%DO_THERMAL_EMISSION
      DO_SURFACE_EMISSION    = self%Surface%Base%TSIO%DO_SURFACE_EMISSION
      DO_D2S_SCALING         = self%Surface%Base%TSIO%DO_D2S_SCALING
!      DO_BRDF_SURFACE 
      DO_USER_OBSGEOMS       = self%Surface%Base%TSIO%DO_USER_OBSGEOMS
      DO_SURFACE_LEAVING     = self%Surface%Base%TSIO%DO_SURFACE_LEAVING
      DO_SL_ISOTROPIC        = self%Surface%Base%TSIO%DO_SL_ISOTROPIC
      DO_PENTADIAG_INVERSE   = self%Surface%Base%TSIO%DO_PENTADIAG_INVERSE

!                         Angles (SZA, viewing, relatuve azimuth), Level
!                        ------------------------------------------------
      HEIGHT_GRID = 0.0
      HEIGHT_GRID(0:NLAYERS) = self.ze * 1.E-3  !HEIGHT_GRID in km
      USER_OBSGEOMS = 0.0
      USER_OBSGEOMS(1,1) = self%Surface%solar_zenith
      USER_OBSGEOMS(1,2) = self%Surface%sensor_zenith
      USER_OBSGEOMS(1,3) = self%Surface%relat_azimuth
      BEAM_SZAS = 0.0
      BEAM_SZAS(1) = self%Surface%solar_zenith
      USER_ANGLES = 0.0
      USER_ANGLES(1) = self%Surface%sensor_zenith
      USER_RELAZMS = 0.0
      USER_RELAZMS(1) = self%Surface%relat_azimuth


!                Populate Scattering Phase Matrix
!                ---------------------------------
! First initialize to zero to be safe
      deltau_vert_input = 0.0
      assym_vert_input  = 0.0
      omega_total_input = 0.0
      d2s_scaling       = 0.0

!                Greek moment for Rayleigh Scattering
!                ----------------------------------------------
      raysmom2 = (1.0 - self%DEPOL_RATIO)/(2.0 + self%DEPOL_RATIO)

!     Loop over the layers:
!     ---------------------
      do i = 1, NLAYERS
         ray_l = self%rot(i)         ! indice l for  each layer
         tau_l = self%tau(i)
         ssa_l = self%ssa(i)
         g_l   = self%g(i)

!        total optical depths for extinction and scattering
!        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         tau_ext = ray_l + tau_l
         tau_scat = ray_l +  ssa_l * tau_l

!        single scattering albedo total
!        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         ssa_tot = tau_scat / tau_ext
         if ( ssa_tot > 0.99999 ) then
            ssa_tot = 0.99999
         end if

         deltau_vert_input(i) = tau_ext
         omega_total_input(i) = ssa_tot

!        assymetry input and d2s_scaling parameter (aerosol + Rayleigh)
!        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!        raysmom(0) = 1
!        raysmom(1) = 0
!        raysmom(2) = raysmom2
!        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!        Add together Aerosol and Rayleigh Parts weighting by scattering optical depth
!        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         aerswt = ssa_l * tau_l / tau_scat
         rayswt = ray_l / tau_scat

         assym_vert_input(i) = g_l*aerswt/tau_scat

         m1 = ray_l * raysmom2
         m2 = 5.0 * g_l * g_l * ssa_l * tau_l
         d2s_scaling(i) = ( m1 + m2 ) /tau_scat / 5.0 

!     end layer loop
!     ---------------
      end do


!     Call the MASTER driver for doing the actual calculation
!     -----------------------------------------------------------
      CALL TWOSTREAM_MASTER &
        ( MAXLAYERS, MAXTOTAL, MAXMESSAGES, MAXBEAMS, MAX_GEOMETRIES,     & 
          MAX_USER_RELAZMS, MAX_USER_ANGLES, MAX_USER_OBSGEOMS,           & 
          DO_UPWELLING, DO_DNWELLING, DO_PLANE_PARALLEL, DO_2S_LEVELOUT,  & 
          DO_MVOUT_ONLY, DO_ADDITIONAL_MVOUT,                             & 
          DO_SOLAR_SOURCES, DO_THERMAL_EMISSION, DO_SURFACE_EMISSION,     & 
          DO_D2S_SCALING, DO_BRDF_SURFACE, DO_USER_OBSGEOMS,              & 
          DO_SURFACE_LEAVING, DO_SL_ISOTROPIC, DO_PENTADIAG_INVERSE,      & 
          self%Surface%Base%TSIO%BVPINDEX,                                &
          self%Surface%Base%TSIO%BVPSCALEFACTOR,                          &
          self%Surface%Base%TSIO%TAYLOR_ORDER,                            &
          TAYLOR_SMALL,                                                   &
          NLAYERS, self%Surface%Base%TSIO%NTOTAL,                         &
          self%Surface%Base%TSIO%STREAM_VALUE,                            &
          self%Surface%Base%N_USER_OBSGEOMS,                              &
          USER_OBSGEOMS,                                                  &
          self%Surface%Base%N_USER_ANGLES,                                &
          USER_ANGLES,                                                    &
          self%Surface%Base%N_USER_RELAZMS,                               &
          USER_RELAZMS,                                                   &
          self%Surface%Base%TSIO%FLUX_FACTOR,                             &
          self%Surface%Base%NBEAMS,                                       &
          BEAM_SZAS,                                                      &
          self%Surface%Base%TSIO%EARTH_RADIUS,                            &
          HEIGHT_GRID,                                                    & 
          DELTAU_VERT_INPUT,                                              & 
          OMEGA_TOTAL_INPUT,                                              &
          ASSYM_VERT_INPUT,                                               &
          D2S_SCALING,                                                    &
          self%Surface%Base%TSIO%THERMAL_BB_INPUT,                        &
          LAMBERTIAN_ALBEDO,                                              &
          self%Surface%Base%TSIO%BRDF_F_0,                                &
          self%Surface%Base%TSIO%BRDF_F,                                  &
          self%Surface%Base%TSIO%UBRDF_F,                                 &
          self%Surface%Base%TSIO%EMISSIVITY,                              &
          self%Surface%Base%TSIO%SURFBB,                                  &
          self%Surface%Base%TSIO%SLTERM_ISOTROPIC,                        &
          self%Surface%Base%TSIO%SLTERM_F_0,                              &
          INTENSITY_TOA, INTENSITY_BOA, FLUXES_TOA, FLUXES_BOA,           & ! Outputs 
          RADLEVEL_UP, RADLEVEL_DN, N_GEOMETRIES,                         & ! Outputs
          self%Surface%Base%TSIO%STATUS_INPUTCHECK,                       & ! Exception handling
          self%Surface%Base%TSIO%C_NMESSAGES,                             & ! Exception handling
          self%Surface%Base%TSIO%C_MESSAGES,                              & ! Exception handling
          self%Surface%Base%TSIO%C_ACTIONS,                               & ! Exception handling
          self%Surface%Base%TSIO%STATUS_EXECUTION,                        & ! Exception handling 
          self%Surface%Base%TSIO%E_MESSAGE,                               & ! Exception handling
          self%Surface%Base%TSIO%E_TRACE_1,                               & ! Exception handling
          self%Surface%Base%TSIO%E_TRACE_2 )                                ! Exception handling

!  Exception handling

        IF ( self%Surface%Base%TSIO%STATUS_INPUTCHECK .eq. 1 ) THEN
          write(*,'(a,i4)')'INPUT Check failed from Baseline Run # ',T
          write(*,*)' - Number of Messages = ', self%Surface%Base%TSIO%C_NMESSAGES
          Do k = 1, self%Surface%Base%TSIO%C_NMESSAGES
            write(*,'(A,I3,A,A)')' - Message # ',K,': ', TRIM(self%Surface%Base%TSIO%C_MESSAGES(K))
            write(*,'(A,I3,A,A)')' - Action  # ',K,': ', TRIM(self%Surface%Base%TSIO%C_ACTIONS(K))
          ENDDO
          stop'Test_2S_only program aborted'
        ENDIF
        IF ( self%Surface%Base%TSIO%STATUS_EXECUTION .eq. 1 ) THEN
          write(*,'(a,i4)')'EXECUTION failed from Baseline Run # ',T
          write(*,*)' - Print 1 Message and 2 Traces'
          write(*,'(A)') TRIM(self%Surface%Base%TSIO%E_MESSAGE)
          write(*,'(A)') TRIM(self%Surface%Base%TSIO%E_TRACE_1)
          write(*,'(A)') TRIM(self%Surface%Base%TSIO%E_TRACE_2)
          stop'Test_2S_only program aborted'
        ENDIF

!     Return TOA radiance
!     -------------------
      output%RADIANCE = INTENSITY_TOA(1)

      output%REFLECTANCE = (pi * output%RADIANCE) / ( cos(self%Surface%solar_zenith/180.0) * self%Surface%Base%TSIO%FLUX_FACTOR)


    end subroutine TWOSTREAM_Run
!.........................................................................
    end module TWOSTREAM_ScatMod
