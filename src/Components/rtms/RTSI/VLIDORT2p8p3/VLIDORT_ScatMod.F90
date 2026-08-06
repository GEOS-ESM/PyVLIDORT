      module VLIDORT_ScatMod
 
      USE VLIDORT_PARS_m
      
      USE VBRDF_SUP_MOD_m
      USE VLIDORT_IO_DEFS_m
      
      USE VLIDORT_AUX_m
      USE VLIDORT_INPUTS_m
      USE VLIDORT_MASTERS_m
      USE vfzmat_Rayleigh_m
      USE vfzmat_Master_m


      USE VLIDORT_Mod
      USE VLIDORT_SurfaceMod
      USE VLIDORT_AOPMod

      implicit NONE

      PUBLIC  VLIDORT_Run      ! Run for each profile (pixel) 
          
      type VLIDORT_output
         real*8, pointer     :: radiance(:)    ! TOA radiance
         real*8, pointer     :: reflectance(:) ! TOA reflectance
         real*8, pointer     :: U(:)           ! U Stokes component
         real*8, pointer     :: Q(:)           ! Q Stokes component
         real*8, pointer     :: V(:)           ! V Stokes component

         real*8, pointer     :: BOA_radiance(:)    ! BOA radiance
         real*8, pointer     :: BOA_reflectance(:) ! BOA reflectance
         real*8, pointer     :: BOA_U(:)           ! U Stokes component
         real*8, pointer     :: BOA_Q(:)           ! Q Stokes component
         real*8, pointer     :: BOA_V(:)           ! V Stokes component

         real*8, pointer     :: ADJUSTED_SLEAVE(:) ! Adjusted water leaving radiance
      end type VLIDORT_output

       
      Contains
!.............................................................................
      subroutine VLIDORT_Run (self, output, rc)
!
!     Computes radiances for a single wavelength, pixel. Optical properties
!     and met fields in self are assumed to have been updated with the
!     apropriate values.
!
      USE VLIDORT_PARS_m
         
      USE VLIDORT_IO_DEFS_m
      USE VBRDF_SUP_MOD_m
      USE VSLEAVE_SUP_MOD_m
      
      USE VLIDORT_AUX_m
      USE VLIDORT_INPUTS_m
      USE VLIDORT_MASTERS_m

      USE VLIDORT_AOPMod

      type(VLIDORT_scat),    intent(inout)        :: self        ! Contains most input
      type(VLIDORT_output),  intent(out)          :: output      ! contains output
      integer,                     intent(out)    :: rc

!                           ----

      integer              :: STATUS_INPUTCHECK, STATUS_CALCULATION 

!                           ----

!     local variables
!     ---------------
      integer                                            :: i, j, k, l, m, n
      integer                                            :: IDR, ierror
      integer                                            :: NLAYERS
 
      logical                                            :: DO_LAMBERTIAN_SURFACE
      real*8                                             :: LAMBERTIAN_ALBEDO

      real*8, dimension(MAX_USER_LEVELS, MAX_GEOMETRIES, &
              MAXSTOKES, MAX_DIRECTIONS)                 :: STOKES
      real*8                                             :: FLUX_FACTOR

      real*8, parameter                                  :: pi = 4.*atan(1.0)
       
      
      rc = 0
 
      if ( .not. self%Surface%Base%initialized ) then
        rc = 1
        return
      end if
      
!                     Stokes/streams/layers/moments
!                     -----------------------------  
      self%Surface%Base%VIO%VLIDORT_FixIn%Cont%TS_NSTOKES                = self%NSTOKES
      self%Surface%Base%VIO%VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT  = self%nmom

      if ( self%nmom .GT. MAXMOMENTS_INPUT)  then
         rc = 5
         return
      end if
      if ( self%NSTOKES  .GT. MAXSTOKES  )   then
         rc = 2 
         return
      end if
      NLAYERS                                                            = self%Surface%Base%VIO%VLIDORT_FixIn%Cont%TS_NLAYERS      
      
      if ( self%Surface%Base%VIO%VLIDORT_FixIn%Bool%TS_DO_UPWELLING ) then
         IDR = 1
      else
         IDR = 2
      end if

!                          Lambertian OR BRDF surface
!                          --------------------------

      if ( self%Surface%sfc_type == LAMBERTIAN ) then               ! Lambertian surface ?
         DO_LAMBERTIAN_SURFACE = .true.        
         LAMBERTIAN_ALBEDO = self%Surface%albedo           ! Lambertian (isotropic) input albedo
      else                                   
         DO_LAMBERTIAN_SURFACE = .false. 
         LAMBERTIAN_ALBEDO = 0.0                           ! Use BRDF -> albedo set to 0
 
      end if

      if ( self%Surface%sfc_type == NOBM_ISO ) then               ! Surface Leaving from External dataset containing isotropic radiance?
         self%Surface%Base%VIO%VLIDORT_Fixin%Bool%TS_DO_SURFACE_LEAVING = .true.
         self%Surface%Base%VIO%VLIDORT_FixIn%Bool%TS_DO_WATER_LEAVING   = .true.
         self%Surface%Base%VIO%VLIDORT_Fixin%Bool%TS_DO_SL_ISOTROPIC    = .true.
         self%Surface%Base%VIO%VLIDORT_Sup%SLEAVE%TS_SLTERM_ISOTROPIC   = self%Surface%Base%VIO%VSLEAVE_Sup_Out%SL_SLTERM_ISOTROPIC
      end if     

      if ( self%Surface%sfc_type == NOBM ) then               ! Surface Leaving from External dataset containing non-isotropic radiance?
         self%Surface%Base%VIO%VLIDORT_Fixin%Bool%TS_DO_SURFACE_LEAVING = .true.
         self%Surface%Base%VIO%VLIDORT_FixIn%Bool%TS_DO_WATER_LEAVING   = .true.
         self%Surface%Base%VIO%VLIDORT_Fixin%Bool%TS_DO_SL_ISOTROPIC    = .false.
         self%Surface%Base%VIO%VLIDORT_Sup%SLEAVE%TS_SLTERM_USERANGLES   = self%Surface%Base%VIO%VSLEAVE_Sup_Out%SL_SLTERM_USERANGLES
      end if 

      if (self%Surface%sleave_adjust) then
          self%Surface%Base%VIO%VLIDORT_Fixin%Bool%TS_DO_TF_ITERATION = .true.
          self%Surface%Base%VIO%VLIDORT_FixIn%Bool%TS_DO_WLADJUSTED_OUTPUT = .true.
      end if

      self%Surface%Base%VIO%VLIDORT_Fixin%Bool%TS_DO_LAMBERTIAN_SURFACE = DO_LAMBERTIAN_SURFACE
      self%Surface%Base%VIO%VLIDORT_FixIn%Optical%TS_LAMBERTIAN_ALBEDO = LAMBERTIAN_ALBEDO
      self%Surface%Base%VIO%VLIDORT_Sup%BRDF%TS_BRDF_F_0        = self%Surface%Base%VIO%VBRDF_Sup_Out%BS_BRDF_F_0
      self%Surface%Base%VIO%VLIDORT_Sup%BRDF%TS_BRDF_F          = self%Surface%Base%VIO%VBRDF_Sup_Out%BS_BRDF_F
      self%Surface%Base%VIO%VLIDORT_Sup%BRDF%TS_USER_BRDF_F_0   = self%Surface%Base%VIO%VBRDF_Sup_Out%BS_USER_BRDF_F_0
      self%Surface%Base%VIO%VLIDORT_Sup%BRDF%TS_USER_BRDF_F     = self%Surface%Base%VIO%VBRDF_Sup_Out%BS_USER_BRDF_F
      self%Surface%Base%VIO%VLIDORT_Sup%BRDF%TS_EXACTDB_BRDFUNC = self%Surface%Base%VIO%VBRDF_Sup_Out%BS_DBOUNCE_BRDFUNC
      self%Surface%Base%VIO%VLIDORT_Sup%BRDF%TS_EMISSIVITY      = self%Surface%Base%VIO%VBRDF_Sup_Out%BS_EMISSIVITY
      self%Surface%Base%VIO%VLIDORT_Sup%BRDF%TS_USER_EMISSIVITY = self%Surface%Base%VIO%VBRDF_Sup_Out%BS_USER_EMISSIVITY

!                         Angles (SZA, viewing, relatuve azimuth), Level
!                        ------------------------------------------------      
      self%Surface%Base%VIO%VLIDORT_ModIn%MSunRays%TS_SZANGLES(1:self%Surface%Base%N_USER_OBSGEOMS) = self%Surface%solar_zenith
      self%Surface%Base%VIO%VLIDORT_ModIn%MUserVal%TS_USER_RELAZMS(1:self%Surface%Base%N_USER_OBSGEOMS) = self%Surface%relat_azimuth
      self%Surface%Base%VIO%VLIDORT_ModIn%MUserVal%TS_USER_VZANGLES_INPUT(1:self%Surface%Base%N_USER_OBSGEOMS) = self%Surface%sensor_zenith
      self%Surface%Base%VIO%VLIDORT_ModIn%MUserVal%TS_N_USER_OBSGEOMS          = self%Surface%Base%N_USER_OBSGEOMS
      self%Surface%Base%VIO%VLIDORT_ModIn%MUserVal%TS_USER_OBSGEOMS_INPUT(1:self%Surface%Base%N_USER_OBSGEOMS,1) = self%Surface%solar_zenith
      self%Surface%Base%VIO%VLIDORT_ModIn%MUserVal%TS_USER_OBSGEOMS_INPUT(1:self%Surface%Base%N_USER_OBSGEOMS,2) = self%Surface%sensor_zenith
      self%Surface%Base%VIO%VLIDORT_ModIn%MUserVal%TS_USER_OBSGEOMS_INPUT(1:self%Surface%Base%N_USER_OBSGEOMS,3) = self%Surface%relat_azimuth      
      self%Surface%Base%VIO%VLIDORT_ModIn%MUserVal%TS_USER_LEVELS(1) = 0.0     ! This is TOA   
      if ( self%DO_BOA ) then
         self%Surface%Base%VIO%VLIDORT_FixIn%Bool%TS_DO_DNWELLING     = .true.
         self%Surface%Base%VIO%VLIDORT_ModIn%MUserVal%TS_USER_LEVELS(2) = 1.0 
      end if

      if (ANY(self%Surface%solar_zenith == 0.0)) then
         self%Surface%Base%VIO%VLIDORT_ModIn%MBool%TS_DO_FOCORR_NADIR        = .true.
         self%Surface%Base%VIO%VLIDORT_ModIn%MBool%TS_DO_FOCORR_OUTGOING     = .false.
      end if           
     
      self%Surface%Base%VIO%VLIDORT_FixIn%Chapman%TS_height_grid(0:NLAYERS)      = self%ze * 1.E-3  ! en km
      self%Surface%Base%VIO%VLIDORT_FixIn%Chapman%TS_pressure_grid(0:NLAYERS)    = self%pe * 1.E-2  ! en hPa
      self%Surface%Base%VIO%VLIDORT_FixIn%Chapman%TS_temperature_grid(0:NLAYERS) = self%te


!                         Get total atmosphere optical inputs
!                        -----------------------------------------------
      call VLIDORT_CombineAOP(self,rc)
      if ( rc /= 0 ) return

      self%Surface%Base%VIO%VLIDORT_FixIn%Optical%TS_DELTAU_VERT_INPUT = self%AOP%deltau_vert_input
      self%Surface%Base%VIO%VLIDORT_ModIn%MOptical%TS_OMEGA_TOTAL_INPUT = self%AOP%omega_total_input
      self%Surface%Base%VIO%VLIDORT_FixIn%Optical%TS_GREEKMAT_TOTAL_INPUT = self%AOP%greekmat_total_input

      if ( self%Surface%Base%VIO%VLIDORT_ModIn%MBool%TS_DO_SSCORR_USEFMAT ) then  ! provide the fmatrices
            self%Surface%Base%VIO%VLIDORT_FixIn%Optical%TS_FMATRIX_UP = self%AOP%fmatrix_up 
            self%Surface%Base%VIO%VLIDORT_FixIn%Optical%TS_FMATRIX_DN = self%AOP%fmatrix_dn
      end if

!     Call the MASTER driver for doing the actual calculation
!     -----------------------------------------------------------
     call VLIDORT_MASTER (self%Surface%Base%DO_DEBUG_INPUT, &
           self%Surface%Base%VIO%VLIDORT_FixIn, &
           self%Surface%Base%VIO%VLIDORT_ModIn, &
           self%Surface%Base%VIO%VLIDORT_Sup, &
           self%Surface%Base%VIO%VLIDORT_Out)
  
      if ( self%Surface%Base%VIO%VLIDORT_Out%Status%TS_STATUS_INPUTCHECK /= 0 ) then
         rc = 3
         write(*,*) 'VLIDORT_MASTER STATUS_INPUTCHECK RETURNED ERROR'
         write(*,*) self%Surface%Base%VIO%VLIDORT_Out%Status%TS_STATUS_INPUTCHECK
         write(*,*) self%Surface%Base%VIO%VLIDORT_Out%Status%TS_CHECKMESSAGES
      end if

      if ( self%Surface%Base%VIO%VLIDORT_Out%Status%TS_STATUS_CALCULATION /= 0 ) then
         write(*,*) 'VLIDORT_MASTER STATUS_CALCULATION RETURNED ERROR'
         write(*,*) self%Surface%Base%VIO%VLIDORT_Out%Status%TS_STATUS_CALCULATION
         write(*,*) self%Surface%Base%VIO%VLIDORT_Out%Status%TS_MESSAGE
         write(*,*) self%Surface%Base%VIO%VLIDORT_Out%Status%TS_TRACE_1
         rc = 4
      end if
      if ( rc /= 0 ) return 

      STOKES = self%Surface%Base%VIO%VLIDORT_Out%Main%TS_STOKES ! output of VLIDORT_MASTER subroutine
      FLUX_FACTOR = self%Surface%Base%VIO%VLIDORT_FixIn%SunRays%TS_FLUX_FACTOR
      
!     Return TOA radiance
!     -------------------
      ! allocate output arrays
      allocate(output%RADIANCE(self%Surface%Base%N_USER_OBSGEOMS))
      allocate(output%REFLECTANCE(self%Surface%Base%N_USER_OBSGEOMS))
      allocate(output%Q(self%Surface%Base%N_USER_OBSGEOMS))
      allocate(output%U(self%Surface%Base%N_USER_OBSGEOMS))
      if (self%NSTOKES == 4) allocate(output%V(self%Surface%Base%N_USER_OBSGEOMS))
      output%RADIANCE    = 0.0
      output%REFLECTANCE = 0.0
      output%Q           = 0
      output%U           = 0
      if (self%NSTOKES == 4) output%V           = 0

      output%RADIANCE = STOKES(1, 1:self%Surface%Base%N_USER_OBSGEOMS, 1, IDR)
      output%Q        = STOKES(1, 1:self%Surface%Base%N_USER_OBSGEOMS, 2, IDR)
      output%U        = STOKES(1, 1:self%Surface%Base%N_USER_OBSGEOMS, 3, IDR)
      if (self%NSTOKES == 4)  output%V = STOKES(1, 1:self%Surface%Base%N_USER_OBSGEOMS, 4, IDR)

      output%REFLECTANCE = (pi * output%RADIANCE) / ( cos(self%Surface%Base%VIO%VLIDORT_ModIn%MSunRays%TS_SZANGLES(1:self%Surface%Base%N_USER_OBSGEOMS)*pi/180.0) * FLUX_FACTOR )   

      if ( self%DO_BOA ) then
         ! allocate output arrays
         allocate(output%BOA_RADIANCE(self%Surface%Base%N_USER_OBSGEOMS))
         allocate(output%BOA_REFLECTANCE(self%Surface%Base%N_USER_OBSGEOMS))
         allocate(output%BOA_Q(self%Surface%Base%N_USER_OBSGEOMS))
         allocate(output%BOA_U(self%Surface%Base%N_USER_OBSGEOMS))
         if (self%NSTOKES == 4) allocate(output%BOA_V(self%Surface%Base%N_USER_OBSGEOMS))
         output%BOA_RADIANCE    = 0.0
         output%BOA_REFLECTANCE = 0.0
         output%BOA_Q           = 0
         output%BOA_U           = 0
         if (self%NSTOKES == 4) output%BOA_V           = 0

         output%BOA_RADIANCE = STOKES(2, 1:self%Surface%Base%N_USER_OBSGEOMS, 1, 2)
         output%Q        = STOKES(2, 1:self%Surface%Base%N_USER_OBSGEOMS, 2, 2)
         output%U        = STOKES(2, 1:self%Surface%Base%N_USER_OBSGEOMS, 3, 2)
         if (self%NSTOKES == 4)  output%BOA_V = STOKES(2, 1:self%Surface%Base%N_USER_OBSGEOMS, 4, 2)

         output%BOA_REFLECTANCE = (pi * output%BOA_RADIANCE) / ( cos(self%Surface%Base%VIO%VLIDORT_ModIn%MSunRays%TS_SZANGLES(1:self%Surface%Base%N_USER_OBSGEOMS)*pi/180.0) * FLUX_FACTOR )   
      end if

      if (self%Surface%sleave_adjust) then
          allocate(output%ADJUSTED_SLEAVE(self%Surface%Base%N_USER_OBSGEOMS))

          output%ADJUSTED_SLEAVE = self%Surface%Base%VIO%VLIDORT_Out%WLOut%TS_WLADJUSTED_ISOTROPIC(1,1:self%Surface%Base%N_USER_OBSGEOMS)
      end if
    
      end subroutine VLIDORT_Run

!.............................................................................

      end module VLIDORT_ScatMod
