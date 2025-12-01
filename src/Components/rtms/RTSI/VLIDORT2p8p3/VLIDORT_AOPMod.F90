      module VLIDORT_AOPMod
 
      USE VLIDORT_PARS_m
      USE vfzmat_Rayleigh_m
      USE vfzmat_Master_m
      
      USE VLIDORT_IO_DEFS_m
      
      USE VLIDORT_AUX_m

      USE VLIDORT_Mod
      USE VLIDORT_SurfaceMod

      implicit NONE

      type VLIDORT_aop
         real*8, pointer          :: greekmat_total_input(:,:,:)
         real*8, pointer          :: deltau_vert_input(:)
         real*8, pointer          :: omega_total_input(:)
         real*8, pointer          :: fmatrix_up(:,:,:)
         real*8, pointer          :: fmatrix_dn(:,:,:)

      end type VLIDORT_aop


      type VLIDORT_scat
         integer         :: NSTOKES             ! Number of stokes vectors
         logical         :: DO_2OS_CORRECTION = .false.   ! Flag to control 2OS Correction (not used here but needed so drivers can work with multiple VLIDORT versions)
         logical         :: DO_BOA = .false.       ! Flag to control whether to do additional down welling calc at BOA
         real*8          :: wavelength          ! in [nm]
         integer         :: nMom                ! number of momemts read (phase function)
         integer         :: nPol                ! number of components of the scattering matrix
         integer         :: N_InAngles          ! number of scattering matrix angles
         real*8, pointer :: InAngles(:)         ! scattering matrix angles
         real*8          :: MISSING             ! MISSING VALUE
         real*8, pointer :: rot(:)              ! rayleigh optical thickness
         real*8, pointer :: depol_ratio
         real*8, pointer :: alpha(:)            ! trace gas absorption optical thickness
         real*8, pointer :: tau(:)              ! aerosol tau
         real*8, pointer :: ssa(:)              ! aerosol ssa
         real*8, pointer ::  pe(:)              ! pressure    at layer edges [Pa]
         real*8, pointer ::  ze(:)              ! height      at layer edges [m]
         real*8, pointer ::  te(:)              ! temperature at layer edges [K]
         real*8, pointer :: pmom(:,:,:)         ! moments of the scattering phase matrix
         real*8, pointer :: fmatrix(:,:,:)      ! scattering phase matrix
         real*8, pointer :: tauI(:)             ! ice cloud tau
         real*8, pointer :: ssaI(:)             ! ice cloud ssa
         real*8, pointer :: pmomI(:,:,:)        ! ice cloud moments of the scattering phase matrix
         real*8, pointer :: fmatrixI(:,:,:)     ! ice cloud scattering phase matrix
         real*8, pointer :: tauL(:)             ! liquid cloud tau
         real*8, pointer :: ssaL(:)             ! liquid cloud ssa
         real*8, pointer :: pmomL(:,:,:)        ! liquid cloud moments of the scattering phase matrix
         real*8, pointer :: fmatrixL(:,:,:)     ! liquid cloud scattering phase matrix

         type(VLIDORT_Surface)           :: Surface
         type(VLIDORT_AOP)               :: AOP

      end type VLIDORT_scat


      PUBLIC  VLIDORT_Rayleigh        ! Calculates layer rayleigh optical thickness
      PUBLIC  VLIDORT_CombineAOP      ! Combines layer atmosphere optical properties to get final VLIDORT inputs
          
       
      Contains

!.............................................................................
      subroutine VLIDORT_Rayleigh (self, rc)
!     Computes Rayleigh optical thickness - populates rot vector in self.  Self contains atmospheric data  
!     Rayleigh extinction profile from Bodhaine et al., (1999) 
!     and Tomasi et al., (2005)
!     (wavelength in micrometer, pressure in hpa)
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    

      USE VLIDORT_PARS_m

      type(VLIDORT_scat), intent(inout)   :: self        ! Contains most input
      integer,            intent(out)     :: rc

      real*8, dimension(0:MAXLAYERS)                     :: Vol    ! Volume coefficient for molecular scattering
      real*8, target, dimension(MAXLAYERS)               :: Ray    ! Layer Rayleigh optical thickness
      real*8, target                                     :: depol

      real*8, dimension(0:MAXLAYERS)                     :: height_grid                     
      real*8, dimension(0:MAXLAYERS)                     :: pressure_grid 
      real*8, dimension(0:MAXLAYERS)                     :: temperature_grid   
      real*8                                             :: wmicron  
      real*8                                             :: Ns
      real*8                                             :: ROD 
      real*8                                             :: ma, C, P_surface
      real*8                                             :: difz
      integer                                            :: NLAYERS
      integer                                            :: j
      real*8, parameter :: Av  = 6.0221367e23    !Avogadros number molecules/mol
      real*8, parameter :: Ts = 288.15  !K, standard temperature
      real*8, parameter :: Ps = 1013.25 !hPa, standard pressture
      real*8, parameter :: g = 980.665  ! gravity in cm/s^2

      NLAYERS                     = self%Surface%Base%VIO%VLIDORT_FixIn%Cont%TS_NLAYERS
      height_grid(0:NLAYERS)      = self%ze * 1.E-3  ! en km
      pressure_grid(0:NLAYERS)    = self%pe * 1.E-2 ! en hPa
      temperature_grid(0:NLAYERS) = self%te         ! K


      wmicron = self%wavelength * 1.E-3  ! micrometer
   ! molecular density for standard pressure and temperature in molecules/cm^3
      Ns = (Av/22.4141)*(273.15/Ts)*(1.0/1000)
      
      do j = 0, NLAYERS 
         Vol(j) = A(wmicron) * Ns * (pressure_grid(j)/Ps) * (Ts/temperature_grid(j))
         Vol(j) = Vol(j) * 1.E5 ! en km-1
      end do
       
!     Simple Interpolation of Volume Scattering Coefficient 
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     
      do j = 0, NLAYERS-1
         difz = height_grid(j) - height_grid(j+1)
         ray(j+1) = difz * (Vol(j) + Vol(j+1))/2.
      end do


      self%rot => Ray
      depol = coef_depol(wmicron)
      self%depol_ratio => depol

      rc = 0

      Contains

!-----------------------------------
! Rayleigh scattering cross section A(wmicron) (cm^2/molecule)
!-----------------------------------      
    
      function A(X)
!  for Rayleigh, X = lambda in micron
          implicit none
          real*8 :: A
          real*8, intent(in) :: X
          real*8 :: Ns, pi3
          real*8 :: XX, XX2, C
          real*8 :: CM
          real*8 :: RI, RI2
          real*8, parameter :: pi = (4.*atan(1.0))
          real*8, parameter :: Av  = 6.0221367e23    !Avogadros number
          real*8, parameter :: Ts = 288.15  !K, standard temperature
         XX=1./X

         XX2=XX*XX

         CM = X/10000.0  ! lambda in cm

         pi3 = pi*pi*pi
        
   !  (ns-1)*1e8  - refractive index of dry air at standard pressure and temperature for CO2 = 300 ppmv
         RI = 8060.51 + 2480990.0/(132.274-XX2) + 17455.7/(39.32957-XX2)
         RI = RI*1e-8 + 1.0

   ! adjust for CO2 = 360 ppm
         C = 360.0*1e-6 ! parts per volume
         RI = (RI-1)*(1.0 + 0.54*(C - 0.0003)) + 1
         RI2 = RI*RI

   ! molecular density for standard pressure and temperature in molecules/cm^3
         Ns = (Av/22.4141)*(273.15/Ts)*(1.0/1000)
            
         A = 24.0*pi3*F_air(XX2,C)* ((RI2-1)*(RI2-1))/(Ns*Ns*CM*CM*CM*CM*(RI2+2)*(RI2+2))
         
         return

      end function A  

      end subroutine VLIDORT_Rayleigh

!.................................................................................

! Rayleigh Depolarization Coefficient
!-----------------------------------

      function coef_depol(X)
         implicit None
         real*8 :: coef_depol
         real*8, intent(in) :: X !lambda in micron
         real*8 :: XX, XX2, C

         XX=1./X

         XX2=XX*XX
   ! CO2 = 360 ppm
         C = 360.0*1e-6 ! parts per volume

         coef_depol = 6.0*(F_air(XX2,C)-1.0)/(7.0*F_air(XX2,C)+3.0)

      end function coef_depol

      function F_air(XX2,C)  
! King factor - depolarization
         implicit None
         real*8 :: F_air
         real*8, intent(in) :: XX2  ! 1/lambda^2 in micron
         real*8, intent(in) :: C    ! CO2 in parts per volume
         real*8 :: depol1, depol2, depol3, depol4

         depol1 = 1.034 + 3.17E-4 * XX2  ! N2
         depol2 = 1.096 + 1.385E-3 * XX2 + 1.448E-4 * XX2 * XX2  ! O2
         depol3 = 1.  !Ar
         depol4 = 1.15  !CO2

         ! C*100  conc in percent parts per volume
         F_air= (78.084*depol1 + 20.946*depol2 + 0.934*depol3 +  C*100.0*depol4) / (78.084 + 20.946 + 0.934 + C*100.0)

      end function F_air       

!.............................................................................
      subroutine VLIDORT_CombineAOP(self, rc)
!
!     Computes radiances for a single wavelength, pixel. Optical properties
!     and met fields in self are assumed to have been updated with the
!     apropriate values.
!
      USE VLIDORT_PARS_m
      USE vfzmat_Rayleigh_m
      USE vfzmat_Master_m
         
      type(VLIDORT_scat),    intent(inout)  :: self        ! Contains most input and container for combined AOPs
      integer,               intent(out)    :: rc


!     local variables
!     ---------------
      integer                                            :: i, j, k, l, m, n
      integer                                            :: IDR, ierror
      integer                                            :: NLAYERS
      real*8                                             :: DEPOL_RATIO
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
 
      real*8, dimension(0:2, 16)                         :: rayvmoms
      real*8, dimension(0:MAXMOMENTS_INPUT,MAXLAYERS,16) :: aervmoms
      real*8, dimension(0:MAXMOMENTS_INPUT,MAXLAYERS,16) :: clLvmoms 
      real*8, dimension(0:MAXMOMENTS_INPUT,MAXLAYERS,16) :: clIvmoms  

      real*8, dimension(MAX_GEOMETRIES,MAXLAYERS,6)      :: rayFmatrices_up
      real*8, dimension(MAX_GEOMETRIES,MAXLAYERS,6)      :: rayFmatrices_dn
      real*8, dimension(0:2, 6)                          :: RayCoeffs
      real*8, dimension(MAX_GEOMETRIES,MAXLAYERS,6)      :: aerFmatrices_up
      real*8, dimension(MAX_GEOMETRIES,MAXLAYERS,6)      :: aerFmatrices_dn
      real*8, dimension(MAXLAYERS, 0:MAXMOMENTS_INPUT, 6):: AerCoeffs
      real*8, dimension(MAX_GEOMETRIES,MAXLAYERS,6)      :: clLFmatrices_up
      real*8, dimension(MAX_GEOMETRIES,MAXLAYERS,6)      :: clLFmatrices_dn
      real*8, dimension(MAXLAYERS, 0:MAXMOMENTS_INPUT, 6):: clLCoeffs
      real*8, dimension(MAX_GEOMETRIES,MAXLAYERS,6)      :: clIFmatrices_up
      real*8, dimension(MAX_GEOMETRIES,MAXLAYERS,6)      :: clIFmatrices_dn
      real*8, dimension(MAXLAYERS, 0:MAXMOMENTS_INPUT, 6):: clICoeffs

      real*8                                             :: difz

      real*8, target, dimension(0:MAXMOMENTS_INPUT,MAXLAYERS,16) :: greekmat_total_input                     
      real*8, target, dimension(MAXLAYERS)                       :: deltau_vert_input
      real*8, target, dimension(MAXLAYERS)                       :: omega_total_input
      real*8, target, dimension(MAXLAYERS, MAX_GEOMETRIES, 6)    :: fmatrix_up
      real*8, target, dimension(MAXLAYERS, MAX_GEOMETRIES, 6)    :: fmatrix_dn

      integer                                            :: N_GEOMS, N_SZAS, N_VZAS, N_AZMS
      integer                                            :: OFFSETS ( MAX_SZANGLES, MAX_USER_VZANGLES )
      real*8                                             :: SZAS (MAX_SZANGLES)
      real*8                                             :: VZAS (MAX_USER_VZANGLES)
      real*8                                             :: AZMS (MAX_USER_RELAZMS)
      real*8                                             :: OBSGEOMS (MAX_GEOMETRIES,3)
      logical                                            :: Exist_InFmatrices ( MAXLAYERS )
      logical                                            :: do_upwelling, do_dnwelling
      logical                                            :: do_ObsGeoms
      logical, parameter                                 :: do_Sunlight = .true.
      real*8                                             :: RayZmatrices_up(MAX_GEOMETRIES,4,4) !not used right now, but required for vfzmat
      real*8                                             :: RayZmatrices_dn(MAX_GEOMETRIES,4,4) !not used right now, but required for vfzmat
      real*8                                             :: AerZmatrices_up(MAX_GEOMETRIES, MAXLAYERS, 4, 4) ! not used right now, but required for vfzmat
      real*8                                             :: AerZmatrices_dn(MAX_GEOMETRIES, MAXLAYERS, 4, 4) ! not used right now, but required for vfzmat
      real*8                                             :: clLZmatrices_up(MAX_GEOMETRIES, MAXLAYERS, 4, 4) ! not used right now, but required for vfzmat
      real*8                                             :: clLZmatrices_dn(MAX_GEOMETRIES, MAXLAYERS, 4, 4) ! not used right now, but required for vfzmat
      real*8                                             :: clIZmatrices_up(MAX_GEOMETRIES, MAXLAYERS, 4, 4) ! not used right now, but required for vfzmat
      real*8                                             :: clIZmatrices_dn(MAX_GEOMETRIES, MAXLAYERS, 4, 4) ! not used right now, but required for vfzmat


      rc = 0
 
      
!                     moments/layers
!                     ---------------  

      if ( self%nmom .GT. MAXMOMENTS_INPUT)  then
         rc = 5
         return
      end if

      NLAYERS  = self%Surface%Base%VIO%VLIDORT_FixIn%Cont%TS_NLAYERS      
      
!                VFZMAT Helper Varibles
!                ----------------------

      do_upwelling = self%Surface%Base%VIO%VLIDORT_FixIn%Bool%TS_DO_UPWELLING
      do_dnwelling = self%Surface%Base%VIO%VLIDORT_FixIn%Bool%TS_DO_DNWELLING
      OFFSETS = 0
      N_GEOMS = self%Surface%Base%VIO%VLIDORT_ModIn%MUserVal%TS_N_USER_OBSGEOMS
      N_SZAS = 0 ; N_VZAS = 0 ; N_AZMS = 0
      SZAS = zero ; VZAS = zero ; AZMS = zero
      Exist_InFmatrices = .false.
      Exist_InFmatrices(1:NLAYERS) = .true.
      do_ObsGeoms = self%Surface%Base%VIO%VLIDORT_ModIn%MBool%TS_DO_OBSERVATION_GEOMETRY   ! set in Mod, must be true
      OBSGEOMS = 0.0 ! initialize for safety
      OBSGEOMS(1:self%Surface%Base%N_USER_OBSGEOMS,:) = self%Surface%Base%VIO%VLIDORT_ModIn%MUserVal%TS_USER_OBSGEOMS_INPUT(1:self%Surface%Base%N_USER_OBSGEOMS,:)

!                Populate Scattering Phase Matrix
!                ---------------------------------


! DEPOL_RATIO is a function of wavelength
      DEPOL_RATIO = self%depol_ratio

      if ( self%Surface%Base%VIO%VLIDORT_ModIn%MBool%TS_DO_SSCORR_USEFMAT ) then  ! Call the vfzmat code
        ! First initialize to zero to be safe
          RayFmatrices_up = 0.0
          RayFmatrices_dn = 0.0
          AerFmatrices_up = 0.0
          AerFmatrices_dn = 0.0
          clLFmatrices_up = 0.0
          clLFmatrices_dn = 0.0
          clIFmatrices_up = 0.0
          clIFmatrices_dn = 0.0
          RayZmatrices_up = 0.0
          RayZmatrices_dn = 0.0
          AerZmatrices_up = 0.0
          AerZmatrices_dn = 0.0
          clLZmatrices_up = 0.0
          clLZmatrices_dn = 0.0
          clIZmatrices_up = 0.0
          clIZmatrices_dn = 0.0
          rayCoeffs       = 0.0
          aerCoeffs       = 0.0
          clLCoeffs       = 0.0
          clICoeffs       = 0.0
          fmatrix_up      = 0.0
          fmatrix_dn      = 0.0
          rayvmoms = 0.0
          aervmoms = 0.0
          clLvmoms = 0.0
          clIvmoms = 0.0

        !  Call to the supplement Rayleigh-master
          Call vfzmat_Rayleigh &
           ( MAX_GEOMETRIES, MAX_SZANGLES, MAX_USER_VZANGLES, MAX_USER_RELAZMS, & ! Input  Dimensions (VLIDORT)
             do_upwelling, do_dnwelling, do_ObsGeoms, do_Sunlight,              & ! Input  Flags
             self%NSTOKES, N_GEOMS, N_SZAS, N_VZAS, N_AZMS,                     & ! Input  Numbers
             OFFSETS, DEG_TO_RAD, SZAS, VZAS, AZMS, OBSGEOMS, DEPOL_RATIO,      & ! Input  Geometries + Depol
             RayFmatrices_up, RayFmatrices_dn, RayZmatrices_up, RayZmatrices_dn, rayCoeffs )

        if (any(self%tau > 0.0)) then
        !  Call to the supplement master for aerosols
          Call vfzmat_Master &
           ( MAXMOMENTS_INPUT, MAX_GEOMETRIES, MAX_SZANGLES, MAX_USER_VZANGLES,        & ! Input  Dimensions (VLIDORT)
             MAX_USER_RELAZMS, MAXLAYERS,                                              & ! Input  Dimensions (VLIDORT)
             self%Surface%Base%Max_InAngles, self%N_InAngles, self%InAngles, self%fmatrix, Exist_InFmatrices,       & ! Input  Fmatrices
             do_upwelling, do_dnwelling, do_ObsGeoms, do_Sunlight,                     & ! Input  Flags
             self%nmom, NLAYERS, self%NSTOKES, N_GEOMS, N_SZAS, N_VZAS, N_AZMS,        & ! Input  Numbers
             OFFSETS, DEG_TO_RAD, SZAS, VZAS, AZMS, OBSGEOMS,                          & ! Input  Geometries
             AerFmatrices_up, AerFmatrices_dn, AerZmatrices_up, AerZmatrices_dn, aerCoeffs )
        end if

        if (any(self%tauL > 0.0)) then
        !  Call to the supplement master for liquid clouds
          Call vfzmat_Master &
           ( MAXMOMENTS_INPUT, MAX_GEOMETRIES, MAX_SZANGLES, MAX_USER_VZANGLES,        & ! Input  Dimensions (VLIDORT)
             MAX_USER_RELAZMS, MAXLAYERS,                                              & ! Input  Dimensions (VLIDORT)
             self%Surface%Base%Max_InAngles, self%N_InAngles, self%InAngles, self%fmatrixL, Exist_InFmatrices,       & ! Input  Fmatrices
             do_upwelling, do_dnwelling, do_ObsGeoms, do_Sunlight,                     & ! Input  Flags
             self%nmom, NLAYERS, self%NSTOKES, N_GEOMS, N_SZAS, N_VZAS, N_AZMS,        & ! Input  Numbers
             OFFSETS, DEG_TO_RAD, SZAS, VZAS, AZMS, OBSGEOMS,                          & ! Input  Geometries
             clLFmatrices_up, clLFmatrices_dn, clLZmatrices_up, clLZmatrices_dn, clLCoeffs )
        end if

        if (any(self%tauI > 0.0)) then
        !  Call to the supplement master for ice clouds
          Call vfzmat_Master &
           ( MAXMOMENTS_INPUT, MAX_GEOMETRIES, MAX_SZANGLES, MAX_USER_VZANGLES,        & ! Input  Dimensions (VLIDORT)
             MAX_USER_RELAZMS, MAXLAYERS,                                              & ! Input  Dimensions (VLIDORT)
             self%Surface%Base%Max_InAngles, self%N_InAngles, self%InAngles, self%fmatrixI, Exist_InFmatrices,       & ! Input  Fmatrices
             do_upwelling, do_dnwelling, do_ObsGeoms, do_Sunlight,                     & ! Input  Flags
             self%nmom, NLAYERS, self%NSTOKES, N_GEOMS, N_SZAS, N_VZAS, N_AZMS,        & ! Input  Numbers
             OFFSETS, DEG_TO_RAD, SZAS, VZAS, AZMS, OBSGEOMS,                          & ! Input  Geometries
             clIFmatrices_up, clIFmatrices_dn, clIZmatrices_up, clIZmatrices_dn, clICoeffs )
        end if



        !  Copy the FMatCoeffs to the right locations of the greekmat
          rayvmoms(0:2,1)  =     rayCoeffs(0:2,1)
          rayvmoms(0:2,2)  = -1.*rayCoeffs(0:2,5)
          rayvmoms(0:2,5)  = -1.*rayCoeffs(0:2,5)
          rayvmoms(0:2,6)  =     rayCoeffs(0:2,2)
          rayvmoms(0:2,11) =     rayCoeffs(0:2,3)
          rayvmoms(0:2,12) = -1.*rayCoeffs(0:2,6)
          rayvmoms(0:2,15) =     rayCoeffs(0:2,6)
          rayvmoms(0:2,16) =     rayCoeffs(0:2,4)

        !     Loop over the layers:
        !     ---------------------
          do i = 1, NLAYERS
        !        Phase function moments - Aerosol and Cloud Part
        !        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
             do l= 0, self%nmom-1
                  aervmoms(l,i,1)  =     aerCoeffs(i,l,1)
                  aervmoms(l,i,2)  = -1.*aerCoeffs(i,l,5)
                  aervmoms(l,i,5)  = -1.*aerCoeffs(i,l,5)
                  aervmoms(l,i,6)  =     aerCoeffs(i,l,2)
                  aervmoms(l,i,11) =     aerCoeffs(i,l,3)
                  aervmoms(l,i,12) = -1.*aerCoeffs(i,l,6)
                  aervmoms(l,i,15) =     aerCoeffs(i,l,6)
                  aervmoms(l,i,16) =     aerCoeffs(i,l,4)           

                  clLvmoms(l,i,1)  =     clLCoeffs(i,l,1)
                  clLvmoms(l,i,2)  = -1.*clLCoeffs(i,l,5)
                  clLvmoms(l,i,5)  = -1.*clLCoeffs(i,l,5)
                  clLvmoms(l,i,6)  =     clLCoeffs(i,l,2)
                  clLvmoms(l,i,11) =     clLCoeffs(i,l,3)
                  clLvmoms(l,i,12) = -1.*clLCoeffs(i,l,6)
                  clLvmoms(l,i,15) =     clLCoeffs(i,l,6)
                  clLvmoms(l,i,16) =     clLCoeffs(i,l,4) 

                  clIvmoms(l,i,1)  =     clICoeffs(i,l,1)
                  clIvmoms(l,i,2)  = -1.*clICoeffs(i,l,5)
                  clIvmoms(l,i,5)  = -1.*clICoeffs(i,l,5)
                  clIvmoms(l,i,6)  =     clICoeffs(i,l,2)
                  clIvmoms(l,i,11) =     clICoeffs(i,l,3)
                  clIvmoms(l,i,12) = -1.*clICoeffs(i,l,6)
                  clIvmoms(l,i,15) =     clICoeffs(i,l,6)
                  clIvmoms(l,i,16) =     clICoeffs(i,l,4)
             end do
          end do

      else  ! use legacy way of constructing VLIDORT scattering coefficients inputs


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
        !                  ----------------------------------------------
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

        !        Aerosol and Cloud phase function moments
        !        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


        !     Loop over the layers:
        !     ---------------------
          do i = 1, NLAYERS
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

          end do
      end if 

!     Combine Gas, Aerosol and Clouds to get final combined atmosphere optical properties
!     ------------------------------------------------------------------------------------------


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
      
!        Add together Aerosol, Cloud, and Rayleigh Parts weighting by scattering optical depth
!        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         aerswt = ssa_l * tau_l / tau_scat  
         rayswt = ray_l / tau_scat 
         clLswt = ssaL_l * tauL_l/ tau_scat
         clIswt = ssaI_l * tauI_l/ tau_scat

         if ( self%Surface%Base%VIO%VLIDORT_ModIn%MBool%TS_DO_SSCORR_USEFMAT ) then  ! Combine the fmatrices
            fmatrix_up(i,1:N_GEOMS,:) = rayswt * RayFmatrices_up(1:N_GEOMS,i,:) + aerswt * AerFmatrices_up(1:N_GEOMS,i,:) + &
                                        clLswt * clLFmatrices_up(1:N_GEOMS,i,:) + clIswt * clIFmatrices_up(1:N_GEOMS,i,:)
            fmatrix_dn(i,1:N_GEOMS,:) = rayswt * RayFmatrices_dn(1:N_GEOMS,i,:) + aerswt * AerFmatrices_dn(1:N_GEOMS,i,:) + &
                                        clLswt * clLFmatrices_dn(1:N_GEOMS,i,:) + clIswt * clIFmatrices_dn(1:N_GEOMS,i,:)
         end if 
 
         ! Both approaches need to combine the moments
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

      self%AOP%deltau_vert_input => deltau_vert_input
      self%AOP%omega_total_input => omega_total_input
      self%AOP%greekmat_total_input => greekmat_total_input
      if ( self%Surface%Base%VIO%VLIDORT_ModIn%MBool%TS_DO_SSCORR_USEFMAT ) then  ! provide the fmatrices
          self%AOP%fmatrix_up => fmatrix_up
          self%AOP%fmatrix_dn => fmatrix_dn
      end if
               
      end subroutine VLIDORT_CombineAOP

!.............................................................................

      end module VLIDORT_AOPMod
