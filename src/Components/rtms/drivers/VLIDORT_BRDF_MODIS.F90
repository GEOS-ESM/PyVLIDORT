module VLIDORT_BRDF_MODIS   
!
!  Simple f77 wrapper for the Python interface to VLIDORT for OMI aerosol
!  channels.
!
!.............................................................................
  implicit NONE

  PUBLIC VLIDORT_LandMODIS_pmom
  PUBLIC VLIDORT_LandMODIS_pmatrix

  contains

  logical function IS_MISSING(x,MISSING)
    real*8, intent(in)     :: x
    real*8, intent(in)     :: MISSING

    IS_MISSING = abs(x/MISSING-1)<0.001
    return
  end function IS_MISSING

  !..........................................................................
  subroutine VLIDORT_LandMODIS_pmom (km, nch, nobs, ngeom, channels, nstreams, do_fullrad, plane_parallel, nMom,  &
                     nPol, NSTOKES, ROT, depol, alpha, tau, ssa, pmom, tauI, ssaI, pmomI, tauL, ssaL, pmomL, &
                     pe, he, te, kernel_wt, param, &
                     solar_zenith, relat_azymuth, sensor_zenith, flux_factor, &
                     MISSING,verbose, debug, radiance_VL_SURF,reflectance_VL_SURF, BR, Q, U, BR_Q, BR_U, rc, &
                     DO_2OS_CORRECTION, DO_BOA)
  !
  ! Place holder.
  !
     use VLIDORT_ScatMod
   
     implicit NONE
 
    logical, parameter                      :: scalar  = .true.   ! is the surface BRDF kernal scalar
    integer, parameter                      :: nkernel = 3   ! number of kernels
    integer, parameter                      :: nparam  = 2   ! number of kernel parameters

  ! !INPUT PARAMETERS:

    integer,          intent(in)            :: km    ! number of vertical levels 
    integer,          intent(in)            :: nch   ! number of channels
    integer,          intent(in)            :: nobs  ! number of observations
    integer,          intent(in)            :: ngeom  ! number of geometries

    logical,          intent(in)            :: plane_parallel ! do plane parallel flag
    logical,          intent(in)            :: do_fullrad ! do both SS and MS, if false does SS only

    integer, target,  intent(in)            :: nMom  ! number of phase function moments 
    integer, target,  intent(in)            :: nPol  ! number of scattering matrix components    
    integer,          intent(in)            :: nstreams  ! number of half space streams
    integer,          intent(in)            :: NSTOKES   ! number of stokes vectors. 1 =scalar 3 or 4 =vector

                    
    real*8, target,   intent(in)            :: channels(nch)    ! wavelengths [nm]

!                                                     ! --- Rayleigh Parameters ---
    real*8, target,   intent(in)            :: ROT(km,nobs,nch) ! rayleigh optical thickness
    real*8, target,   intent(in)            :: depol(nch)       ! rayleigh depolarization ratio used in phase matrix

!                                                   ! --- Trace Gas Absorption ---
    real*8, target,   intent(in)            :: alpha(km,nch,nobs) ! trace gas absoprtion optical thickness

  !                                                   ! --- Aerosol Optical Properties ---
    real*8, target,   intent(in)            :: tau(km,nch,nobs) ! aerosol optical depth
    real*8, target,   intent(in)            :: ssa(km,nch,nobs) ! single scattering albedo    
    real*8, target,   intent(in)            :: pmom(km,nch,nobs,nMom,nPol) !components of the scat phase matrix

    real*8, target,   intent(in)            :: tauI(km,nch,nobs) ! ice cloud optical depth
    real*8, target,   intent(in)            :: ssaI(km,nch,nobs) ! ice cloud single scattering albedo   
    real*8, target,   intent(in)            :: pmomI(km,nch,nobs,nMom,nPol) !ice cloud components of the scat phase matrix

    real*8, target,   intent(in)            :: tauL(km,nch,nobs) ! liquid cloud optical depth
    real*8, target,   intent(in)            :: ssaL(km,nch,nobs) ! liquid cloud single scattering albedo   
    real*8, target,   intent(in)            :: pmomL(km,nch,nobs,nMom,nPol) !liquid cloud components of the scat phase matrix

    real*8, target,   intent(in)            :: MISSING          ! MISSING VALUE
    real*8, target,   intent(in)            :: pe(km+1,nobs)    ! pressure at layer edges [Pa]
    real*8, target,   intent(in)            :: he(km+1,nobs)    ! height above sea-level  [m]
    real*8, target,   intent(in)            :: te(km+1,nobs)    ! temperature at layer edges [K]

    real*8, target,   intent(in)            :: kernel_wt(nkernel,nch,nobs)    ! kernel weights (/fiso,fgeo,fvol/)
    real*8, target,   intent(in)            :: param(nparam,nch,nobs)         ! Li-Sparse parameters 
                                                                              ! param1 = crown relative height (h/b)
                                                                              ! param2 = shape parameter (b/r)
                         
    real*8, target,   intent(in)            :: solar_zenith(nobs,ngeom)  
    real*8, target,   intent(in)            :: relat_azymuth(nobs,ngeom) 
    real*8, target,   intent(in)            :: sensor_zenith(nobs,ngeom) 

    real*8,           intent(in)            :: flux_factor(nch) ! solar flux (F0)
    integer,          intent(in)            :: verbose
    logical,          intent(in)            :: debug

  ! !OUTPUT PARAMETERS:
    real*8,           intent(out)           :: radiance_VL_SURF(nobs,nch,ngeom)     ! TOA normalized radiance from VLIDORT using surface module
    real*8,           intent(out)           :: reflectance_VL_SURF(nobs, nch, ngeom) ! TOA reflectance from VLIDORT using surface module
    integer,          intent(out)           :: rc                             ! return code

    real*8,           intent(out)           :: BR(nobs,nch,ngeom)                   ! bidirectional reflectance 
    real*8,           intent(out)           :: BR_Q(nobs,nch,ngeom)                 ! bidirectional reflectance Q 
    real*8,           intent(out)           :: BR_U(nobs,nch,ngeom)                 ! bidirectional reflectance U
    real*8,           intent(out)           :: Q(nobs, nch, ngeom)                   ! Stokes parameter Q
    real*8,           intent(out)           :: U(nobs, nch, ngeom)                   ! Stokes parameter U   

  ! !OPTIONAL PARAMETERS
    logical, optional, intent(in) :: DO_2OS_CORRECTION
    logical, optional, intent(in) :: DO_BOA

  !                               ---
    
    integer             :: i,j,n,p,ier,g

    
    type(VLIDORT_scat) :: SCAT
    type(VLIDORT_output)  :: output

    rc = 0
    ier = 0
   
    if (present(DO_BOA)) SCAT%DO_BOA = DO_BOA
    SCAT%Surface%Base%NSTREAMS = nstreams
    SCAT%Surface%Base%DO_FULLRAD_MODE   = do_fullrad
    SCAT%Surface%Base%DO_PLANE_PARALLEL = plane_parallel
    SCAT%Surface%Base%DO_DEBUG_INPUT    = debug
    SCAT%Surface%Base%N_USER_OBSGEOMS   = ngeom
    SCAT%Surface%Base%NBEAMS            = ngeom
    SCAT%Surface%Base%N_USER_STREAMS    = ngeom
    SCAT%Surface%Base%N_USER_RELAZMS    = ngeom
    SCAT%Surface%Base%NGREEK_MOMENTS_INPUT = nMom
    SCAT%Surface%Base%USEFMAT = .false.
    call VLIDORT_Init( SCAT%Surface%Base, km, rc, SCAT%DO_BOA)
    if ( rc /= 0 ) return

    SCAT%nMom    = nMom
    SCAT%nPol    = nPol
    if (present(DO_2OS_CORRECTION)) then
      SCAT%DO_2OS_CORRECTION = DO_2OS_CORRECTION
      if (DO_2OS_CORRECTION) then
        SCAT%NSTOKES = 1
      else
        SCAT%NSTOKES = NSTOKES
      end if
    else    
      SCAT%NSTOKES = NSTOKES
    end if

    if ( SCAT%NSTOKES  .GT. MAXSTOKES  )   return

    do j = 1, nobs
       
      ! Make sure angles are available
      ! -----------------------------------------
      if ( IS_MISSING(solar_zenith(j,1),MISSING)  .OR. & 
           IS_MISSING(sensor_zenith(j,1),MISSING) .OR. &
           IS_MISSING(relat_azymuth(j,1),MISSING)  )  then

        radiance_VL_SURF(j,:,:) = MISSING
        reflectance_VL_SURF(j,:,:) = MISSING
        Q(j,:,:) = MISSING
        U(j,:,:) = MISSING
        BR(j,:,:) = MISSING
        BR_Q(j,:,:) = MISSING
        BR_U(j,:,:) = MISSING 
        cycle

      end if

      SCAT%pe => pe(:,j)
      SCAT%ze => he(:,j)
      SCAT%te => te(:,j) 

      do i = 1, nch
        ! set solar flux
        SCAT%Surface%Base%VIO%VLIDORT_FixIn%SunRays%TS_FLUX_FACTOR = flux_factor(i)
      
        ! Make sure kernel weights and parameters are defined
        do n = 1, nkernel
          if (IS_MISSING(kernel_wt(n,i,j),MISSING)) then
            radiance_VL_SURF(j,i,:) = MISSING
            reflectance_VL_SURF(j,i,:) = MISSING
            Q(j,i,:) = MISSING
            U(j,i,:) = MISSING
            BR(j,i,:) = MISSING
            BR_Q(j,i,:) = MISSING
            BR_U(j,i,:) = MISSING
            cycle
          end if
        end do
        do p = 1, nparam
          if (IS_MISSING(param(p,i,j),MISSING)) then
            radiance_VL_SURF(j,i,:) = MISSING
            reflectance_VL_SURF(j,i,:) = MISSING
            Q(j,i,:) = MISSING
            U(j,i,:) = MISSING
            BR(j,i,:) = MISSING
            BR_Q(j,i,:) = MISSING
            BR_U(j,i,:) = MISSING
            cycle
          end if
        end do

        if ( verbose > 0 ) then
          print*, 'DO MODIS BRDF'
        end if

        call VLIDORT_LANDMODIS(SCAT%Surface,solar_zenith(j,:),&
                               sensor_zenith(j,:),relat_azymuth(j,:),&
                               kernel_wt(1,i,j),kernel_wt(2,i,j),kernel_wt(3,i,j),&
                               reshape(param(:,i,j),(/nparam/)),&
                               scalar,rc)

        SCAT%wavelength = channels(i) 
        SCAT%rot => rot(:,j,i)  
        SCAT%depol_ratio => depol(i)    
        SCAT%alpha => alpha(:,i,j) 
        SCAT%tau => tau(:,i,j)
        SCAT%ssa => ssa(:,i,j)
        SCAT%pmom => pmom(:,i,j,:,:)
        SCAT%tauI => tauI(:,i,j)
        SCAT%ssaI => ssaI(:,i,j)
        SCAT%pmomI => pmomI(:,i,j,:,:)    
        SCAT%tauL => tauL(:,i,j)
        SCAT%ssaL => ssaL(:,i,j)
        SCAT%pmomL => pmomL(:,i,j,:,:)            

        do g = 1,ngeom
          BR(j,i,g) = SCAT%Surface%Base%VIO%VBRDF_Sup_Out%BS_DBOUNCE_BRDFUNC(1,1,1,g)
          BR_Q(j,i,g) = 0
          BR_U(j,i,g) = 0     

          ! Check to make sure kernel combination is not unphysical
          if (BR(j,i,g) < 0) then
            radiance_VL_SURF(j,i,g)    = -500
            reflectance_VL_SURF(j,i,g) = -500
            BR(j,i,g)    = -500
            BR_Q(j,i,g)    = -500
            BR_U(j,i,g)    = -500
            Q(j,i,g)     = -500
            U(j,i,g)     = -500
            cycle
          end if
        end do
        call VLIDORT_Run (SCAT, output, ier)

        if (SCAT%DO_BOA) then
            radiance_VL_SURF(j,i,:)    = output%BOA_radiance
            reflectance_VL_SURF(j,i,:) = output%BOA_reflectance
            if (NSTOKES .eq. 1) then
                ! this should be done upstream,  but impose here explicitly just in case
                Q(j,i,:)                   = 0.0
                U(j,i,:)                   = 0.0
            else
                Q(j,i,:)                   = output%BOA_Q
                U(j,i,:)                   = output%BOA_U               
            end if 
        else
            radiance_VL_SURF(j,i,:)    = output%radiance
            reflectance_VL_SURF(j,i,:) = output%reflectance
            if (NSTOKES .eq. 1) then
                ! this should be done upstream,  but impose here explicitly just in case
                Q(j,i,:)                   = 0.0
                U(j,i,:)                   = 0.0
            else
                Q(j,i,:)                   = output%Q
                U(j,i,:)                   = output%U                
            end if
        end if

        if ( ier /= 0 ) then
            radiance_VL_SURF(j,i,:) = MISSING
            reflectance_VL_SURF(j,i,:) = MISSING               
            Q(j,i,:) = MISSING
            U(j,i,:) = MISSING
            BR(j,i,:) = MISSING
            BR_Q(j,i,:) = MISSING
            BR_U(j,i,:) = MISSING
            cycle
        end if

      end do ! end loop over channels
       
      if ( verbose > 0 ) then
        if ( mod(j-1,1000) == 0 ) then
          print *, '<> VLIDORT Vector: ', nint(j*100./nobs), '%'
        end if
      end if

    end do ! Loop over obs

  end subroutine VLIDORT_LandMODIS_pmom


  !..........................................................................
  subroutine VLIDORT_LandMODIS_pmatrix (km, nch, nobs, ngeom, channels, nstreams, do_fullrad, plane_parallel, nAng,  &
                     nPol, NSTOKES, InAngles, ROT, depol, alpha, tau, ssa, pmatrix, tauI, ssaI, pmatrixI, tauL, ssaL, pmatrixL, &
                     pe, he, te, kernel_wt, param, &
                     solar_zenith, relat_azymuth, sensor_zenith, flux_factor, &
                     MISSING,verbose, debug, radiance_VL_SURF,reflectance_VL_SURF, BR, Q, U, BR_Q, BR_U, rc, &
                     DO_2OS_CORRECTION, DO_BOA)
  !
  ! Place holder.
  !
     use VLIDORT_ScatMod
   
     implicit NONE
 
    logical, parameter                      :: scalar  = .true.   ! is the surface BRDF kernal scalar
    integer, parameter                      :: nkernel = 3   ! number of kernels
    integer, parameter                      :: nparam  = 2   ! number of kernel parameters

  ! !INPUT PARAMETERS:

    integer,          intent(in)            :: km    ! number of vertical levels 
    integer,          intent(in)            :: nch   ! number of channels
    integer,          intent(in)            :: nobs  ! number of observations
    integer,          intent(in)            :: ngeom  ! number of geometries

    logical,          intent(in)            :: plane_parallel ! do plane parallel flag
    logical,          intent(in)            :: do_fullrad ! do both SS and MS, if false do SS only

    integer, target,  intent(in)            :: nAng  ! number of phase function angles
    integer, target,  intent(in)            :: nPol  ! number of scattering matrix components    
    integer,          intent(in)            :: nstreams  ! number of half space streams
    integer,          intent(in)            :: NSTOKES   ! number of stokes vectors. 1 =scalar 3 or 4 =vector

                    
    real*8, target,   intent(in)            :: channels(nch)    ! wavelengths [nm]
    real*8, target,   intent(in)            :: InAngles(nAng)   ! phase matrix angles

!                                                     ! --- Rayleigh Parameters ---
    real*8, target,   intent(in)            :: ROT(km,nobs,nch) ! rayleigh optical thickness
    real*8, target,   intent(in)            :: depol(nch)       ! rayleigh depolarization ratio used in phase matrix

!                                                   ! --- Trace Gas Absorption ---
    real*8, target,   intent(in)            :: alpha(km,nch,nobs) ! trace gas absoprtion optical thickness

  !                                                   ! --- Aerosol Optical Properties ---
    real*8, target,   intent(in)            :: tau(km,nch,nobs) ! aerosol optical depth
    real*8, target,   intent(in)            :: ssa(km,nch,nobs) ! single scattering albedo    
    real*8, target,   intent(in)            :: pmatrix(km,nch,nobs,nAng,nPol) !components of the scat phase matrix

    real*8, target,   intent(in)            :: tauI(km,nch,nobs) ! ice cloud optical depth
    real*8, target,   intent(in)            :: ssaI(km,nch,nobs) ! ice cloud single scattering albedo   
    real*8, target,   intent(in)            :: pmatrixI(km,nch,nobs,nAng,nPol) !ice cloud components of the scat phase matrix

    real*8, target,   intent(in)            :: tauL(km,nch,nobs) ! liquid cloud optical depth
    real*8, target,   intent(in)            :: ssaL(km,nch,nobs) ! liquid cloud single scattering albedo   
    real*8, target,   intent(in)            :: pmatrixL(km,nch,nobs,nAng,nPol) !liquid cloud components of the scat phase matrix

    real*8, target,   intent(in)            :: MISSING          ! MISSING VALUE
    real*8, target,   intent(in)            :: pe(km+1,nobs)    ! pressure at layer edges [Pa]
    real*8, target,   intent(in)            :: he(km+1,nobs)    ! height above sea-level  [m]
    real*8, target,   intent(in)            :: te(km+1,nobs)    ! temperature at layer edges [K]

    real*8, target,   intent(in)            :: kernel_wt(nkernel,nch,nobs)    ! kernel weights (/fiso,fgeo,fvol/)
    real*8, target,   intent(in)            :: param(nparam,nch,nobs)         ! Li-Sparse parameters 
                                                                              ! param1 = crown relative height (h/b)
                                                                              ! param2 = shape parameter (b/r)
                         
    real*8, target,   intent(in)            :: solar_zenith(nobs,ngeom)  
    real*8, target,   intent(in)            :: relat_azymuth(nobs,ngeom) 
    real*8, target,   intent(in)            :: sensor_zenith(nobs,ngeom) 

    real*8,           intent(in)            :: flux_factor(nch) ! solar flux (F0)
    integer,          intent(in)            :: verbose
    logical,          intent(in)            :: debug

  ! !OUTPUT PARAMETERS:
    real*8,           intent(out)           :: radiance_VL_SURF(nobs,nch,ngeom)     ! TOA normalized radiance from VLIDORT using surface module
    real*8,           intent(out)           :: reflectance_VL_SURF(nobs, nch, ngeom) ! TOA reflectance from VLIDORT using surface module
    integer,          intent(out)           :: rc                             ! return code

    real*8,           intent(out)           :: BR(nobs,nch,ngeom)                   ! bidirectional reflectance 
    real*8,           intent(out)           :: BR_Q(nobs,nch,ngeom)                 ! bidirectional reflectance Q 
    real*8,           intent(out)           :: BR_U(nobs,nch,ngeom)                 ! bidirectional reflectance U
    real*8,           intent(out)           :: Q(nobs, nch, ngeom)                   ! Stokes parameter Q
    real*8,           intent(out)           :: U(nobs, nch, ngeom)                   ! Stokes parameter U   

  ! !OPTIONAL PARAMETERS
    logical, optional, intent(in) :: DO_2OS_CORRECTION
    logical, optional, intent(in) :: DO_BOA

  !                               ---
    
    integer             :: i,j,n,p,ier,g

    
    type(VLIDORT_scat) :: SCAT
    type(VLIDORT_output)  :: output

    rc = 0
    ier = 0
   
    if (present(DO_BOA)) SCAT%DO_BOA = DO_BOA
    SCAT%Surface%Base%NSTREAMS = nstreams
    SCAT%Surface%Base%DO_PLANE_PARALLEL = plane_parallel
    SCAT%Surface%Base%DO_FULLRAD_MODE   = do_fullrad
    SCAT%Surface%Base%DO_DEBUG_INPUT    = debug
    SCAT%Surface%Base%N_USER_OBSGEOMS   = ngeom
    SCAT%Surface%Base%NBEAMS            = ngeom
    SCAT%Surface%Base%N_USER_STREAMS    = ngeom
    SCAT%Surface%Base%N_USER_RELAZMS    = ngeom
    SCAT%Surface%Base%Max_InAngles      = nAng
    SCAT%Surface%Base%USEFMAT = .true.
    call VLIDORT_Init( SCAT%Surface%Base, km, rc, SCAT%DO_BOA)
    if ( rc /= 0 ) return

    SCAT%N_InAngles    = nAng
    SCAT%InAngles => InAngles
    SCAT%nMom    = SCAT%Surface%Base%NGREEK_MOMENTS_INPUT  ! using the default value
    SCAT%nPol    = nPol
    if (present(DO_2OS_CORRECTION)) then
      SCAT%DO_2OS_CORRECTION = DO_2OS_CORRECTION
      if (DO_2OS_CORRECTION) then
        SCAT%NSTOKES = 1
      else
        SCAT%NSTOKES = NSTOKES
      end if
    else    
      SCAT%NSTOKES = NSTOKES
    end if

    if ( SCAT%NSTOKES  .GT. MAXSTOKES  )   return

    do j = 1, nobs
       
      ! Make sure angles are available
      ! -----------------------------------------
      if ( IS_MISSING(solar_zenith(j,1),MISSING)  .OR. & 
           IS_MISSING(sensor_zenith(j,1),MISSING) .OR. &
           IS_MISSING(relat_azymuth(j,1),MISSING)  )  then

        radiance_VL_SURF(j,:,:) = MISSING
        reflectance_VL_SURF(j,:,:) = MISSING
        Q(j,:,:) = MISSING
        U(j,:,:) = MISSING
        BR(j,:,:) = MISSING
        BR_Q(j,:,:) = MISSING
        BR_U(j,:,:) = MISSING 
        cycle

      end if

      SCAT%pe => pe(:,j)
      SCAT%ze => he(:,j)
      SCAT%te => te(:,j) 

      do i = 1, nch
        ! set solar flux
        SCAT%Surface%Base%VIO%VLIDORT_FixIn%SunRays%TS_FLUX_FACTOR = flux_factor(i)
      
        ! Make sure kernel weights and parameters are defined
        do n = 1, nkernel
          if (IS_MISSING(kernel_wt(n,i,j),MISSING)) then
            radiance_VL_SURF(j,i,:) = MISSING
            reflectance_VL_SURF(j,i,:) = MISSING
            Q(j,i,:) = MISSING
            U(j,i,:) = MISSING
            BR(j,i,:) = MISSING
            BR_Q(j,i,:) = MISSING
            BR_U(j,i,:) = MISSING
            cycle
          end if
        end do
        do p = 1, nparam
          if (IS_MISSING(param(p,i,j),MISSING)) then
            radiance_VL_SURF(j,i,:) = MISSING
            reflectance_VL_SURF(j,i,:) = MISSING
            Q(j,i,:) = MISSING
            U(j,i,:) = MISSING
            BR(j,i,:) = MISSING
            BR_Q(j,i,:) = MISSING
            BR_U(j,i,:) = MISSING
            cycle
          end if
        end do

        if ( verbose > 0 ) then
          print*, 'DO MODIS BRDF'
        end if

        call VLIDORT_LANDMODIS(SCAT%Surface,solar_zenith(j,:),&
                               sensor_zenith(j,:),relat_azymuth(j,:),&
                               kernel_wt(1,i,j),kernel_wt(2,i,j),kernel_wt(3,i,j),&
                               reshape(param(:,i,j),(/nparam/)),&
                               scalar,rc)

        SCAT%wavelength = channels(i) 
        SCAT%rot => rot(:,j,i)  
        SCAT%depol_ratio => depol(i)    
        SCAT%alpha => alpha(:,i,j) 
        SCAT%tau => tau(:,i,j)
        SCAT%ssa => ssa(:,i,j)
        SCAT%fmatrix => pmatrix(:,i,j,:,:)
        SCAT%tauI => tauI(:,i,j)
        SCAT%ssaI => ssaI(:,i,j)
        SCAT%fmatrixI => pmatrixI(:,i,j,:,:)    
        SCAT%tauL => tauL(:,i,j)
        SCAT%ssaL => ssaL(:,i,j)
        SCAT%fmatrixL => pmatrixL(:,i,j,:,:)            

        do g = 1,ngeom
          BR(j,i,g) = SCAT%Surface%Base%VIO%VBRDF_Sup_Out%BS_DBOUNCE_BRDFUNC(1,1,1,g)
          BR_Q(j,i,g) = 0
          BR_U(j,i,g) = 0     

          ! Check to make sure kernel combination is not unphysical
          if (BR(j,i,g) < 0) then
            radiance_VL_SURF(j,i,g)    = -500
            reflectance_VL_SURF(j,i,g) = -500
            BR(j,i,g)    = -500
            BR_Q(j,i,g)    = -500
            BR_U(j,i,g)    = -500
            Q(j,i,g)     = -500
            U(j,i,g)     = -500
            cycle
          end if
        end do
        call VLIDORT_Run (SCAT, output, ier)

        if (SCAT%DO_BOA) then
            radiance_VL_SURF(j,i,:)    = output%BOA_radiance
            reflectance_VL_SURF(j,i,:) = output%BOA_reflectance
            if (NSTOKES .eq. 1) then
                ! this should be done upstream,  but impose here explicitly just in case
                Q(j,i,:)                   = 0.0
                U(j,i,:)                   = 0.0
            else
                Q(j,i,:)                   = output%BOA_Q
                U(j,i,:)                   = output%BOA_U               
            end if 
        else
            radiance_VL_SURF(j,i,:)    = output%radiance
            reflectance_VL_SURF(j,i,:) = output%reflectance
            if (NSTOKES .eq. 1) then
                ! this should be done upstream,  but impose here explicitly just in case
                Q(j,i,:)                   = 0.0
                U(j,i,:)                   = 0.0
            else
                Q(j,i,:)                   = output%Q
                U(j,i,:)                   = output%U                
            end if
        end if

        if ( ier /= 0 ) then
            radiance_VL_SURF(j,i,:) = MISSING
            reflectance_VL_SURF(j,i,:) = MISSING               
            Q(j,i,:) = MISSING
            U(j,i,:) = MISSING
            BR(j,i,:) = MISSING
            BR_Q(j,i,:) = MISSING
            BR_U(j,i,:) = MISSING
            cycle
        end if

      end do ! end loop over channels
       
      if ( verbose > 0 ) then
        if ( mod(j-1,1000) == 0 ) then
          print *, '<> VLIDORT Vector: ', nint(j*100./nobs), '%'
        end if
      end if

    end do ! Loop over obs

  end subroutine VLIDORT_LandMODIS_pmatrix



end module VLIDORT_BRDF_MODIS
