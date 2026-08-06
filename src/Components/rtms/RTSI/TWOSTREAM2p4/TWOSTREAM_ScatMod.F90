      module TWOSTREAM_ScatMod

      USE TWOSTREAM_PARS
      USE vfzmat_Rayleigh_m
      USE vfzmat_Master_m
      USE twostream_brdf_supplement_m
      USE twostream_master_m

      USE TWOSTREAM_Mod
      USE TWOSTREAM_SurfaceMod

      implicit NONE

      PUBLIC  TWOSTREAM_Run

      type TWOSTREAM_scat

         real*8          :: wavelength          ! in [nm]
         integer         :: nMom                ! number of momemts read (phase function)
         integer         :: N_InAngles          ! number of scattering matrix angles
         real*8, pointer :: InAngles(:)         ! scattering matrix angles
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
         real*8, pointer :: fmatrix(:,:,:)      ! scattering phase matrix
         real*8, pointer :: tauI(:)             ! ice cloud tau
         real*8, pointer :: ssaI(:)             ! ice cloud ssa
         real*8, pointer :: gI(:)               ! ice cloud asymmetry factor
         real*8, pointer :: fmatrixI(:,:,:)     ! ice cloud scattering phase matrix
         real*8, pointer :: tauL(:)             ! liquid cloud tau
         real*8, pointer :: ssaL(:)             ! liquid cloud ssa
         real*8, pointer :: gL(:)               ! liquid cloud asymmetry factor
         real*8, pointer :: fmatrixL(:,:,:)     ! liquid cloud scattering phase matrix

         type(TWOSTREAM_Surface) :: Surface

      end type TWOSTREAM_scat

      type TWOSTREAM_output
         real*8, pointer     :: RADIANCE(:)    ! TOA radiance
         real*8, pointer     :: REFLECTANCE(:) ! TOA reflectance
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
      USE vfzmat_Rayleigh_m
      USE vfzmat_Master_m
      USE TWOSTREAM_MOD
      USE twostream_brdf_supplement_m

      USE twostream_master_m

      USE twostream_Master_debug_m, only: write_twostream_Master_debug
      USE twostream_Master_outputs_debug_m, only: write_twostream_Master_outputs_debug
      USE vfzmat_Master_debug_m, only: write_vfzmat_Master_debug

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
      real*8                                             :: alpha_l
      real*8                                             :: tauI_l
      real*8                                             :: ssaI_l
      real*8                                             :: gI_l
      real*8                                             :: tauL_l
      real*8                                             :: ssaL_l
      real*8                                             :: gL_l
      real*8                                             :: tau_ext
      real*8                                             :: tau_scat
      real*8                                             :: ssa_tot
      real*8                                             :: raysmom2
      real*8                                             :: aersmom2
      real*8                                             :: clLsmom2
      real*8                                             :: clIsmom2
      real*8                                             :: aerswt
      real*8                                             :: rayswt
      real*8                                             :: clLswt
      real*8                                             :: clIswt

      INTEGER, PARAMETER :: MAXMOMENTS_INPUT = 2
      INTEGER, PARAMETER :: NSTOKES = 1

      real*8, target, dimension(0:MAXMOMENTS_INPUT,MAXLAYERS,16) :: aervmoms
      real*8, dimension(0:MAXMOMENTS_INPUT,MAXLAYERS,16) :: clLvmoms
      real*8, dimension(0:MAXMOMENTS_INPUT,MAXLAYERS,16) :: clIvmoms

      real*8, dimension(MAX_GEOMETRIES,MAXLAYERS,6)      :: aerFmatrices_up
      real*8, dimension(MAX_GEOMETRIES,MAXLAYERS,6)      :: aerFmatrices_dn
      real*8, dimension(MAXLAYERS, 0:MAXMOMENTS_INPUT, 6):: AerCoeffs
      real*8, dimension(MAX_GEOMETRIES,MAXLAYERS,6)      :: clLFmatrices_up
      real*8, dimension(MAX_GEOMETRIES,MAXLAYERS,6)      :: clLFmatrices_dn
      real*8, dimension(MAXLAYERS, 0:MAXMOMENTS_INPUT, 6):: clLCoeffs
      real*8, dimension(MAX_GEOMETRIES,MAXLAYERS,6)      :: clIFmatrices_up
      real*8, dimension(MAX_GEOMETRIES,MAXLAYERS,6)      :: clIFmatrices_dn
      real*8, dimension(MAXLAYERS, 0:MAXMOMENTS_INPUT, 6):: clICoeffs

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

      integer                                            :: N_GEOMS, N_SZAS, N_VZAS, N_AZMS
      integer                                            :: OFFSETS ( MAX_SZANGLES, MAX_USER_ANGLES )
      real*8                                             :: SZAS (MAX_SZANGLES)
      real*8                                             :: VZAS (MAX_USER_ANGLES)
      real*8                                             :: AZMS (MAX_USER_RELAZMS)
      real*8                                             :: OBSGEOMS (MAX_GEOMETRIES,3)
      logical                                            :: Exist_InFmatrices ( MAXLAYERS )
      logical                                            :: do_ObsGeoms
      logical, parameter                                 :: do_Sunlight = .true.
      real*8                                             :: AerZmatrices_up(MAX_GEOMETRIES, MAXLAYERS, 4, 4) ! not used right now, but required for vfzmat
      real*8                                             :: AerZmatrices_dn(MAX_GEOMETRIES, MAXLAYERS, 4, 4) ! not used right now, but required for vfzmat
      real*8                                             :: clLZmatrices_up(MAX_GEOMETRIES, MAXLAYERS, 4, 4) ! not used right now, but required for vfzmat
      real*8                                             :: clLZmatrices_dn(MAX_GEOMETRIES, MAXLAYERS, 4, 4) ! not used right now, but required for vfzmat
      real*8                                             :: clIZmatrices_up(MAX_GEOMETRIES, MAXLAYERS, 4, 4) ! not used right now, but required for vfzmat
      real*8                                             :: clIZmatrices_dn(MAX_GEOMETRIES, MAXLAYERS, 4, 4) ! not used right now, but required for vfzmat

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
      USER_OBSGEOMS(1:self%Surface%Base%N_USER_OBSGEOMS,1) = self%Surface%solar_zenith
      USER_OBSGEOMS(1:self%Surface%Base%N_USER_OBSGEOMS,2) = self%Surface%sensor_zenith
      USER_OBSGEOMS(1:self%Surface%Base%N_USER_OBSGEOMS,3) = self%Surface%relat_azimuth
      BEAM_SZAS = 0.0
      BEAM_SZAS(1:self%Surface%Base%N_USER_OBSGEOMS) = self%Surface%solar_zenith
      USER_ANGLES = 0.0
      USER_ANGLES(1:self%Surface%Base%N_USER_OBSGEOMS) = self%Surface%sensor_zenith
      USER_RELAZMS = 0.0
      USER_RELAZMS(1:self%Surface%Base%N_USER_OBSGEOMS) = self%Surface%relat_azimuth

!                VFZMAT Helper Varibles
!                ----------------------

      OFFSETS = 0
      N_GEOMS = self%Surface%Base%N_USER_OBSGEOMS
      N_SZAS = 1 ; N_VZAS = 1 ; N_AZMS = 1
      Exist_InFmatrices = .false.
      Exist_InFmatrices(1:NLAYERS) = .true.
      do_ObsGeoms = self%Surface%Base%TSIO%DO_USER_OBSGEOMS   ! set in Mod, must be true
      OBSGEOMS = 0.0 ! initialize for safety
      OBSGEOMS(1:self%Surface%Base%N_USER_OBSGEOMS,:) = USER_OBSGEOMS(1:self%Surface%Base%N_USER_OBSGEOMS,:)
      SZAS = OBSGEOMS(1,1)
      VZAS = OBSGEOMS(1,2)
      AZMS = OBSGEOMS(1,3)

!                Populate Scattering Phase Matrix
!                ---------------------------------

! Call the vfzmat code
        ! First initialize to zero to be safe
          AerFmatrices_up = 0.0
          AerFmatrices_dn = 0.0
          clLFmatrices_up = 0.0
          clLFmatrices_dn = 0.0
          clIFmatrices_up = 0.0
          clIFmatrices_dn = 0.0
          AerZmatrices_up = 0.0
          AerZmatrices_dn = 0.0
          clLZmatrices_up = 0.0
          clLZmatrices_dn = 0.0
          clIZmatrices_up = 0.0
          clIZmatrices_dn = 0.0
          aerCoeffs       = 0.0
          clLCoeffs       = 0.0
          clICoeffs       = 0.0

        if (any(self%tau > 0.0)) then
        !  Call to the supplement master for aerosols

          Call vfzmat_Master &
           ( MAXMOMENTS_INPUT, MAX_GEOMETRIES, MAX_SZANGLES, MAX_USER_ANGLES,        & ! Input  Dimensions (VLIDORT)
             MAX_USER_RELAZMS, MAXLAYERS,                                              & ! Input  Dimensions (VLIDORT)
             self%Surface%Base%Max_InAngles, self%N_InAngles, self%InAngles, self%fmatrix, Exist_InFmatrices,       & ! Input  Fmatrices
             do_upwelling, do_dnwelling, do_ObsGeoms, do_Sunlight,                     & ! Input  Flags
             self%nmom, NLAYERS, NSTOKES, N_GEOMS, N_SZAS, N_VZAS, N_AZMS,        & ! Input  Numbers
             OFFSETS, DEG_TO_RAD, SZAS, VZAS, AZMS, OBSGEOMS,                          & ! Input  Geometries
             AerFmatrices_up, AerFmatrices_dn, AerZmatrices_up, AerZmatrices_dn, aerCoeffs )

           if ( self%Surface%Base%DO_DEBUG_INPUT ) then
                ! --- Call the debug writer ---

                call write_vfzmat_Master_debug( &
                    MAXMOMENTS_INPUT, MAX_GEOMETRIES, MAX_SZANGLES, MAX_USER_ANGLES, &
                    MAX_USER_RELAZMS, MAXLAYERS,                                        &
                    self%Surface%Base%Max_InAngles, self%N_InAngles,                   &
                    self%InAngles, self%fmatrix, Exist_InFmatrices,                    &
                    do_upwelling, do_dnwelling, do_ObsGeoms, do_Sunlight,              &
                    self%nmom, NLAYERS, 1, N_GEOMS, N_SZAS, N_VZAS, N_AZMS,&
                    OFFSETS, DEG_TO_RAD, SZAS, VZAS, AZMS, OBSGEOMS,                  &
                    AerFmatrices_up, AerFmatrices_dn,                                  &
                    AerZmatrices_up, AerZmatrices_dn,                                  &
                    AerCoeffs )
           end if
        end if

        if (any(self%tauL > 0.0)) then
        !  Call to the supplement master for liquid clouds
          Call vfzmat_Master &
           ( MAXMOMENTS_INPUT, MAX_GEOMETRIES, MAX_SZANGLES, MAX_USER_ANGLES,        & ! Input  Dimensions (VLIDORT)
             MAX_USER_RELAZMS, MAXLAYERS,                                              & ! Input  Dimensions (VLIDORT)
             self%Surface%Base%Max_InAngles, self%N_InAngles, self%InAngles, self%fmatrixL, Exist_InFmatrices,       & ! Input  Fmatrices
             do_upwelling, do_dnwelling, do_ObsGeoms, do_Sunlight,                     & ! Input  Flags
             self%nmom, NLAYERS, NSTOKES, N_GEOMS, N_SZAS, N_VZAS, N_AZMS,        & ! Input  Numbers
             OFFSETS, DEG_TO_RAD, SZAS, VZAS, AZMS, OBSGEOMS,                          & ! Input  Geometries
             clLFmatrices_up, clLFmatrices_dn, clLZmatrices_up, clLZmatrices_dn, clLCoeffs )
        end if

        if (any(self%tauI > 0.0)) then
        !  Call to the supplement master for ice clouds
          Call vfzmat_Master &
           ( MAXMOMENTS_INPUT, MAX_GEOMETRIES, MAX_SZANGLES, MAX_USER_ANGLES,        & ! Input  Dimensions (VLIDORT)
             MAX_USER_RELAZMS, MAXLAYERS,                                              & ! Input  Dimensions (VLIDORT)
             self%Surface%Base%Max_InAngles, self%N_InAngles, self%InAngles, self%fmatrixI, Exist_InFmatrices,       & ! Input  Fmatrices
             do_upwelling, do_dnwelling, do_ObsGeoms, do_Sunlight,                     & ! Input  Flags
             self%nmom, NLAYERS, NSTOKES, N_GEOMS, N_SZAS, N_VZAS, N_AZMS,        & ! Input  Numbers
             OFFSETS, DEG_TO_RAD, SZAS, VZAS, AZMS, OBSGEOMS,                          & ! Input  Geometries
             clIFmatrices_up, clIFmatrices_dn, clIZmatrices_up, clIZmatrices_dn, clICoeffs )
        end if

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
         alpha_l = self%alpha(i)
         tau_l = self%tau(i)
         ssa_l = self%ssa(i)
         g_l   = self%g(i)
         tauI_l = self%tauI(i)
         tauL_l = self%tauL(i)
         ssaI_l = self%ssaI(i)
         ssaL_l = self%ssaL(i)
         gL_l   = self%gL(i)
         gI_l   = self%gI(i)

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

!        Greek moment for Aerosol and Cloud Scattering
!        ----------------------------------------------
         aersmom2 = aerCoeffs(i,2,1)
         clLsmom2 = clLCoeffs(i,2,1)
         clIsmom2 = clICoeffs(i,2,1)

!        assymetry input and d2s_scaling parameter (aerosol + Rayleigh)
!        d2s_scaling is the second phase function moment
!        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!        raysmom(0) = 1
!        raysmom(1) = 0
!        raysmom(2) = raysmom2
!        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!        Add together Aerosol and Rayleigh Parts weighting by scattering optical depth
!        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         aerswt = ssa_l * tau_l / tau_scat
         rayswt = ray_l / tau_scat
         clLswt = ssaL_l * tauL_l/ tau_scat
         clIswt = ssaI_l * tauI_l/ tau_scat

         assym_vert_input(i) = g_l*aerswt + gL_l*clLswt + gI_l*clIswt

         d2s_scaling(i) = (rayswt*raysmom2 + aerswt*aersmom2 + clLswt*clLsmom2 + clIswt*clIsmom2)/5.0

!     end layer loop
!     ---------------
      end do


!     Call the MASTER driver for doing the actual calculation
!     -----------------------------------------------------------
      if ( self%Surface%Base%DO_DEBUG_INPUT ) then
                ! --- Call the debug writer ---
          CALL write_twostream_master_debug &
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
              self%Surface%Base%TSIO%SLTERM_F_0)
      end if

      CALL TWOSTREAM_MASTER &
        ( MAXLAYERS, MAXTOTAL, MAXMESSAGES, MAXBEAMS, MAX_GEOMETRIES,     & 
          MAX_USER_RELAZMS, MAX_USER_ANGLES, MAX_USER_OBSGEOMS,           & 
          DO_UPWELLING, DO_DNWELLING, DO_PLANE_PARALLEL, DO_2S_LEVELOUT,  & 
          DO_MVOUT_ONLY, DO_ADDITIONAL_MVOUT,                             & 
          DO_SOLAR_SOURCES, DO_THERMAL_EMISSION, DO_SURFACE_EMISSION,     & 
          DO_D2S_SCALING, DO_BRDF_SURFACE, DO_USER_OBSGEOMS,              & 
          DO_SURFACE_LEAVING, DO_SL_ISOTROPIC, DO_PENTADIAG_INVERSE,      &
          self%Surface%Base%DO_DEBUG_INPUT,                               & 
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

      if ( self%Surface%Base%DO_DEBUG_INPUT ) then
          CALL write_twostream_master_outputs_debug &
            ( MAXLAYERS, MAXMESSAGES, MAXBEAMS, MAX_GEOMETRIES,       & 
              NLAYERS, self%Surface%Base%NBEAMS, N_GEOMETRIES,        & 
              INTENSITY_TOA, INTENSITY_BOA, FLUXES_TOA, FLUXES_BOA,   & 
              RADLEVEL_UP, RADLEVEL_DN,                               & 
              self%Surface%Base%TSIO%STATUS_INPUTCHECK, self%Surface%Base%TSIO%C_NMESSAGES, &
              self%Surface%Base%TSIO%C_MESSAGES, self%Surface%Base%TSIO%C_ACTIONS,          &
              self%Surface%Base%TSIO%STATUS_EXECUTION, self%Surface%Base%TSIO%E_MESSAGE,    &
              self%Surface%Base%TSIO%E_TRACE_1, self%Surface%Base%TSIO%E_TRACE_2 )
      end if

!  Exception handling

        IF ( self%Surface%Base%TSIO%STATUS_INPUTCHECK .eq. 1 ) THEN
          write(*,'(a)')'INPUT Check failed'
          write(*,*)' - Number of Messages = ', self%Surface%Base%TSIO%C_NMESSAGES
          Do k = 1, self%Surface%Base%TSIO%C_NMESSAGES
            write(*,'(A,I3,A,A)')' - Message # ',K,': ', TRIM(self%Surface%Base%TSIO%C_MESSAGES(K))
            write(*,'(A,I3,A,A)')' - Action  # ',K,': ', TRIM(self%Surface%Base%TSIO%C_ACTIONS(K))
          ENDDO
          stop'Test_2S_only program aborted'
        ENDIF
        IF ( self%Surface%Base%TSIO%STATUS_EXECUTION .eq. 1 ) THEN
          write(*,'(a)')'EXECUTION failed'
          write(*,*)' - Print 1 Message and 2 Traces'
          write(*,'(A)') TRIM(self%Surface%Base%TSIO%E_MESSAGE)
          write(*,'(A)') TRIM(self%Surface%Base%TSIO%E_TRACE_1)
          write(*,'(A)') TRIM(self%Surface%Base%TSIO%E_TRACE_2)
          stop'Test_2S_only program aborted'
        ENDIF

!     Return TOA radiance
!     -------------------
      allocate(output%RADIANCE(self%Surface%Base%N_USER_OBSGEOMS))
      allocate(output%REFLECTANCE(self%Surface%Base%N_USER_OBSGEOMS))
      output%RADIANCE = INTENSITY_TOA(1:self%Surface%Base%N_USER_OBSGEOMS)

      output%REFLECTANCE = (pi * output%RADIANCE) / ( cos(self%Surface%solar_zenith/180.0) * self%Surface%Base%TSIO%FLUX_FACTOR)


    end subroutine TWOSTREAM_Run
!.........................................................................
    end module TWOSTREAM_ScatMod
