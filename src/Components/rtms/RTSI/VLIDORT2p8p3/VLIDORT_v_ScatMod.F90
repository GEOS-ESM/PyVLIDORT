      module VLIDORT_ScatMod
 
      USE VLIDORT_PARS_m
      
      USE VBRDF_SUP_MOD_m
      USE VLIDORT_IO_DEFS_m
      
      USE VLIDORT_AUX_m
      USE VLIDORT_INPUTS_m
      USE VLIDORT_MASTERS_m

      USE VLIDORT_Mod
      USE VLIDORT_SurfaceMod
      USE VLIDORT_AOPMod

      implicit NONE

      PUBLIC  VLIDORT_Run_Vector      ! Run for each profile (pixel) 
          

      type VLIDORT_scat
         integer         :: NSTOKES             ! Number of stokes vectors
         logical         :: DO_2OS_CORRECTION = .false.   ! Flag to control 2OS Correction (not used here but needed so drivers can work with multiple VLIDORT versions)
         logical         :: DO_BOA = .false.       ! Flag to control whether to do additional down welling calc at BOA
         real*8          :: wavelength          ! in [nm]
         integer         :: nMom                ! number of momemts read (phase function)
         integer         :: nPol                ! number of components of the scattering matrix
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
         real*8, pointer :: pmom(:,:,:)          ! components of the scattering phase matrix
         real*8, pointer :: tauI(:)             ! ice cloud tau
         real*8, pointer :: ssaI(:)             ! ice cloud ssa
         real*8, pointer ::   gI(:)             ! ice cloud asymmetry factor
         real*8, pointer :: pmomI(:,:,:)        ! ice cloud components of the scattering phase matrix
         real*8, pointer :: tauL(:)             ! liquid cloud tau
         real*8, pointer :: ssaL(:)             ! liquid cloud ssa
         real*8, pointer ::   gL(:)             ! liquid cloud asymmetry factor
         real*8, pointer :: pmomL(:,:,:)        ! liquid cloud components of the scattering phase matrix

         type(VLIDORT_Surface) :: Surface
         type(VLIDORT_AOP)               :: AOP

      end type VLIDORT_scat



      type VLIDORT_output_vector
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
      end type VLIDORT_output_vector

       
      Contains
!.............................................................................
      subroutine VLIDORT_Run_Vector (self, output, rc)
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

      type(VLIDORT_scat),    intent(inout)  :: self        ! Contains most input
      type(VLIDORT_output_vector), intent(out)    :: output      ! contains output
      integer,                     intent(out)    :: rc

!                           ----

      integer              :: STATUS_INPUTCHECK, STATUS_CALCULATION 

!                           ----

!     local variables
!     ---------------
      integer                                            :: i, j, k, l, m, n
      integer                                            :: IDR, ierror
      integer                                            :: NLAYERS
      real*8                                             :: ray_l     
      real*8                                             :: alpha_l 
      real*8                                             :: tau_l 
      real*8                                             :: ssa_l       
      real*8                                             :: ssaL_l 
      real*8                                             :: tauL_l 
      real*8                                             :: ssaI_l 
      real*8                                             :: tauI_l    
      real*8                                             :: tau_ext
      real*8                                             :: tau_scat
      real*8                                             :: ssa_tot
      real*8                                             :: raysmom2
      real*8                                             :: gammamom2
      real*8                                             :: alphamom2 
      real*8                                             :: deltamom1 
      real*8                                             :: aerswt 
      real*8                                             :: rayswt
      real*8                                             :: clIswt
      real*8                                             :: clLswt
 
      real*8, dimension(0:MAXMOMENTS_INPUT,MAXLAYERS,16) :: aervmoms
      real*8, dimension(0:MAXMOMENTS_INPUT,MAXLAYERS,16) :: clLvmoms 
      real*8, dimension(0:MAXMOMENTS_INPUT,MAXLAYERS,16) :: clIvmoms  
      real*8, dimension(0:2, 16)                         :: rayvmoms
      real*8                                             :: difz

      logical                                            :: DO_LAMBERTIAN_SURFACE
      real*8                                             :: LAMBERTIAN_ALBEDO

      real*8, dimension(0:MAXMOMENTS_INPUT,MAXLAYERS,16) :: greekmat_total_input                     
      real*8, dimension(MAXLAYERS)                       :: deltau_vert_input
      real*8, dimension(MAXLAYERS)                       :: omega_total_input

      real*8, dimension(MAX_USER_LEVELS, MAX_GEOMETRIES, &
              MAXSTOKES, MAX_DIRECTIONS)                 :: STOKES
      real*8                                             :: FLUX_FACTOR

      real*8, parameter                                  :: pi = 4.*atan(1.0)
      real*8                                             :: DEPOL_RATIO
       
      
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


! DEPOL_RATIO is a function of wavelength in microns
      DEPOL_RATIO = coef_depol(self%wavelength * 1.E-3)
!                Populate Scattering Phase Matrix
!                ---------------------------------
! First initialize to zero to be safe
      rayvmoms = 0.0
      aervmoms = 0.0 
      clLvmoms = 0.0
      clIvmoms = 0.0     

!                Greek moments for Rayleigh Scattering 
!                The same for all layers because DEPOL_RATIO is
!                taken as constant right now. 
!                ----------------------------------------------
      gammamom2 = -SQRT(6.) * (1 - DEPOL_RATIO) / (2 + DEPOL_RATIO)
      alphamom2 = 6 * (1 - DEPOL_RATIO) / (2 + DEPOL_RATIO)
      deltamom1 = 3 * (1 - 2 * DEPOL_RATIO) / (2 + DEPOL_RATIO)
      raysmom2 = (1.0 - DEPOL_RATIO)/(2.0 + DEPOL_RATIO) 
   
      rayvmoms(0,1) = 1.0
      rayvmoms(1,1) = 0.0
      rayvmoms(2,1) = raysmom2
      rayvmoms(0,2) = 0.0
      rayvmoms(1,2) = 0.0
      rayvmoms(2,2) = gammamom2
      do k = 3, 4
         do l = 0, 2
            rayvmoms(l,k) = 0.0
         end do
      end do
      rayvmoms(0,5) = 0.0
      rayvmoms(1,5) = 0.0
      rayvmoms(2,5) = gammamom2 
      rayvmoms(0,6) = 0.0
      rayvmoms(1,6) = 0.0
      rayvmoms(2,6) = alphamom2 
      do k = 7, 15
         do l = 0, 2
            rayvmoms(l,k) = 0.0
         end do
      end do
      rayvmoms(0,16) = 0.0
      rayvmoms(1,16) = deltamom1
      rayvmoms(2,16) = 0.0 

!     Loop over the layers:
!     ---------------------
      do i = 1, NLAYERS  
         ray_l = self%rot(i)         ! indice l for  each layer 
         ! check to see if alpha is initialized - alpha isn't
         ! implemented everywhere yet
         if (associated(self%alpha)) then
             alpha_l = self%alpha(i)
         else
             alpha_l = 0.0
         end if 
         tau_l = self%tau(i)
         ssa_l = self%ssa(i) 
         tauI_l = self%tauI(i)
         tauL_l = self%tauL(i)
         ssaI_l = self%ssaI(i)
         ssaL_l = self%ssaL(i)         
        
!        total optical depths for extinction and scattering 
!        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         tau_ext = alpha_l + ray_l + tau_l + tauL_l + tauI_l
         tau_scat = ray_l +  ssa_l * tau_l + tauL_l*ssaL_l + tauI_l*ssaI_l

!        single scattering albedo total
!        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         ssa_tot = tau_scat / tau_ext
         if ( ssa_tot > 0.99999 ) then
            ssa_tot = 0.99999
         end if
     
         deltau_vert_input(i) = tau_ext
         omega_total_input(i) = ssa_tot 

!        VECTOR phase function moments 
!        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~      
      
!        Phase function moments - Aerosol Part
!        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         do l= 0, self%nmom-1               
            aervmoms(l,i,1)  = self%pmom(i,l+1,1) ! P11 
            aervmoms(l,i,2)  = self%pmom(i,l+1,2) ! P12                
            aervmoms(l,i,5)  = self%pmom(i,l+1,2) ! P12 = P21                            
            aervmoms(l,i,11) = self%pmom(i,l+1,3) ! P33 
            aervmoms(l,i,12) = self%pmom(i,l+1,4) ! P34            
            aervmoms(l,i,15) = -self%pmom(i,l+1,4) ! - P34
            aervmoms(l,i,6)  = self%pmom(i,l+1,5) ! P22     
            aervmoms(l,i,16) = self%pmom(i,l+1,6) ! P44  
         end do

!        Phase function moments - Cloud Part
!        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         do l= 0, self%nmom-1               
            clLvmoms(l,i,1)  = self%pmomL(i,l+1,1) ! P11 
            clIvmoms(l,i,1)  = self%pmomI(i,l+1,1) ! P11 
         end do
         
!        Add together Aerosol, Cloud, and Rayleigh Parts weighting by scattering optical depth
!        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         aerswt = ssa_l * tau_l / tau_scat  
         rayswt = ray_l / tau_scat 
         clLswt = ssaL_l * tauL_l/ tau_scat
         clIswt = ssaI_l * tauI_l/ tau_scat
         do k = 1, 16
            do l = 0,2
               greekmat_total_input(l,i,k) = rayvmoms(l,k) * rayswt + aervmoms(l,i,k) * aerswt + &
                                             clLvmoms(l,i,k) * clLswt + clIvmoms(l,i,k) * clIswt
            end do
        
            do l = 3, self%nmom
               greekmat_total_input(l,i,k) = aervmoms(l,i,k) * aerswt + clLvmoms(l,i,k) * clLswt + clIvmoms(l,i,k) * clIswt
            end do
               
         end do
         greekmat_total_input(0,i,1) = 1.0
      
           
!     end layer loop
!     ---------------
      end do
  
      self%Surface%Base%VIO%VLIDORT_FixIn%Optical%TS_DELTAU_VERT_INPUT = deltau_vert_input
      self%Surface%Base%VIO%VLIDORT_ModIn%MOptical%TS_OMEGA_TOTAL_INPUT = omega_total_input
      self%Surface%Base%VIO%VLIDORT_FixIn%Optical%TS_GREEKMAT_TOTAL_INPUT = greekmat_total_input
       
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
    
      end subroutine VLIDORT_Run_Vector

!.............................................................................

      end module VLIDORT_v_ScatMod
