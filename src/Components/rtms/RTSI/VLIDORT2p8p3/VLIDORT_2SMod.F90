      module VLIDORT_2SMod
 
      USE VLIDORT_PARS_m
      
      USE VLIDORT_IO_DEFS_m
      
      USE VLIDORT_INPUTS_m
      USE vfzmat_Master_m


      USE VLIDORT_Mod
      USE VLIDORT_SurfaceMod
      USE VLIDORT_AOPMod

      implicit NONE

      PUBLIC  VLIDORT_2SInputs      ! Get the TWO STREAM inputs
          
      type VLIDORT_AERSMOM2
         real*8, pointer     :: aersmom2(:,:)    ! aerosol first two moments of expanded scattering phase function
      end type VLIDORT_AERSMOM2

       
      Contains
!.............................................................................
      subroutine VLIDORT_2SInputs (self, output, rc)
!
!     Computed the first two aerosol phase function moments.
!     To be passed on to TWOSTREAM code 
!
      USE VLIDORT_PARS_m
         
      USE VLIDORT_IO_DEFS_m
      
      USE VLIDORT_INPUTS_m

      USE VLIDORT_AOPMod

      type(VLIDORT_scat),    intent(inout)        :: self        ! Contains most input
      type(VLIDORT_AERSMOM2),  intent(out)        :: output      ! contains output
      integer,                     intent(out)    :: rc

!                           ----

      integer              :: STATUS_INPUTCHECK, STATUS_CALCULATION 

!                           ----

!     local variables
!     ---------------
      integer                                            :: i, j, k, l, m, n
      integer                                            :: IDR, ierror
      integer                                            :: NLAYERS
 
      rc = 0
 
      if ( .not. self%Surface%Base%initialized ) then
        rc = 1
        return
      end if
      
!                     Stokes/streams/layers/moments
!                     -----------------------------  
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

      output%aersmom2 = self%AOP%aervmoms(0:1,:,1)
      if ( rc /= 0 ) return

      end subroutine VLIDORT_2SInputs

!.............................................................................

      end module VLIDORT_2SMod
