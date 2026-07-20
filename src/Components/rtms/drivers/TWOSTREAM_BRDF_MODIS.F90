module TWOSTREAM_BRDF_MODIS   
!
!  Simple f77 wrapper for the Python interface to 2STREAM
!
!.............................................................................
  implicit NONE

  PUBLIC TWOSTREAM_BRDF_LandMODIS

  contains

  logical function IS_MISSING(x,MISSING)
    real*8, intent(in)     :: x
    real*8, intent(in)     :: MISSING

    IS_MISSING = abs(x/MISSING-1)<0.001
    return
  end function IS_MISSING

  subroutine TWOSTREAM_BRDF_LandMODIS (km, nch, nobs, ngeom, channels, plane_parallel, nAng, &
                     nPol, InAngles, ROT, depol_ratio, alpha, tau, ssa, g, pmatrix, &
                     tauI, ssaI, gI, pmatrixI, tauL, ssaL, gL, pmatrixL, &
                     pe, he, te, kernel_wt, param, &
                     solar_zenith, relat_azymuth, sensor_zenith, flux_factor, &
                     MISSING,verbose,debug,radiance_L_SURF,reflectance_L_SURF, rc )
  !
  ! Uses TWOSTREAM to compute TOA radiances.
  !
    use TWOSTREAM_ScatMod


    implicit NONE

    integer, parameter            :: nkernel = 3
    integer, parameter            :: nparam  = 2
  ! !INPUT PARAMETERS:
    integer,          intent(in)  :: km    ! number of vertical levels
    integer,          intent(in)  :: nch   ! number of channels
    integer,          intent(in)  :: nobs  ! number of observations
    integer,          intent(in)  :: ngeom  ! number of geometries

    logical,          intent(in)  :: plane_parallel ! do plane parallel flag

    integer, target,  intent(in)  :: nAng  ! number of phase function angles
    integer, target,  intent(in)  :: nPol  ! number of components
                                        
    real*8, target,   intent(in)  :: channels(nch)    ! wavelengths [nm]
    real*8, target,   intent(in)  :: InAngles(nAng)   ! phase matrix angles

  !                                                   ! --- Aerosol Optical Properties ---
    real*8, target,   intent(in)  :: tau(km,nch,nobs) ! aerosol optical depth
    real*8, target,   intent(in)  :: ssa(km,nch,nobs) ! single scattering albedo
    real*8, target,   intent(in)  :: g(km,nch,nobs)   ! asymmetry factor
    real*8, target,   intent(in)  :: pmatrix(km,nch,nobs,nAng,nPol) !components of the scat phase matrix

    real*8, target,   intent(in)  :: tauI(km,nch,nobs) ! ice cloud optical depth
    real*8, target,   intent(in)  :: ssaI(km,nch,nobs) ! ice cloud single scattering albedo
    real*8, target,   intent(in)  :: gI(km,nch,nobs)   ! ice cloud asymmetry factor
    real*8, target,   intent(in)  :: pmatrixI(km,nch,nobs,nAng,nPol) !ice cloud components of the scat phase matrix

    real*8, target,   intent(in)  :: tauL(km,nch,nobs) ! liquid cloud optical depth
    real*8, target,   intent(in)  :: ssaL(km,nch,nobs) ! liquid cloud single scattering albedo
    real*8, target,   intent(in)  :: gL(km,nch,nobs)   ! liquid cloud asymmetry factor
    real*8, target,   intent(in)  :: pmatrixL(km,nch,nobs,nAng,nPol) !liquid cloud components of the scat phase matrix

    real*8, target,   intent(in)  :: alpha(km,nch,nobs) ! trace gas absorption
    real*8, target,   intent(in)  :: ROT(km,nobs,nch) ! rayleigh optical thickness
    real*8, target,   intent(in)  :: depol_ratio(nch) ! depolarization ratio

    real*8, target,   intent(in)  :: pe(km+1,nobs)    ! pressure at layer edges [Pa]
    real*8, target,   intent(in)  :: he(km+1,nobs)    ! height above sea-level  [m]
    real*8, target,   intent(in)  :: te(km+1,nobs)    ! temperature at layer edges [K]

    real*8, target,   intent(in)  :: kernel_wt(nkernel,nch,nobs)   ! kernel weights (/fiso,fgeo,fvol/)
    real*8, target,   intent(in)  :: param(nparam,nch,nobs)        ! Li-Sparse parameters 
                                                                   ! param1 = crown relative height (h/b)
                                                                   ! param2 = shape parameter (b/r)
    
    real*8, target,   intent(in)  :: MISSING          ! MISSING VALUE                                      
    real*8, target,   intent(in)  :: solar_zenith(nobs,ngeom)  
    real*8, target,   intent(in)  :: relat_azymuth(nobs,ngeom) 
    real*8, target,   intent(in)  :: sensor_zenith(nobs,ngeom) 

    real*8,           intent(in)  :: flux_factor(nch) ! solar flux (F0)
    
    integer,          intent(in)  :: verbose
    logical,          intent(in)  :: debug

  ! !OUTPUT PARAMETERS:

    real*8,           intent(out) :: radiance_L_SURF(nobs,nch,ngeom)       ! TOA normalized radiance from LIDORT using surface module
    real*8,           intent(out) :: reflectance_L_SURF(nobs,nch,ngeom)   ! TOA reflectance from LIDORT using surface module
    integer,          intent(out) :: rc                               ! return code

    integer             :: i,j,n,p,ier

  
    type(TWOSTREAM_scat)           :: SCAT
    type(TWOSTREAM_output)         :: output

    rc = 0
    ier = 0

    SCAT%Surface%Base%DO_PLANE_PARALLEL = plane_parallel
    SCAT%Surface%Base%DO_DEBUG_INPUT    = debug
    SCAT%Surface%Base%N_USER_OBSGEOMS   = ngeom
    SCAT%Surface%Base%NBEAMS            = ngeom
    SCAT%Surface%Base%N_USER_ANGLES     = ngeom
    SCAT%Surface%Base%N_USER_RELAZMS    = ngeom

    SCAT%N_InAngles    = nAng
    SCAT%InAngles => InAngles

    call TWOSTREAM_Init( SCAT%Surface%Base, km, rc)
    if ( rc /= 0 ) then
      write(*,*) 'TWOSTREAM_Init returning with error'
      return
    endif 


    SCAT%nMom    = 2

    do j = 1,nobs

       ! Make sure albedo and angles are available
       ! -----------------------------------------
      if ( IS_MISSING(solar_zenith(j,1),MISSING)  .OR. & 
           IS_MISSING(sensor_zenith(j,1),MISSING) .OR. &
           IS_MISSING(relat_azymuth(j,1),MISSING)  )  then

        radiance_L_SURF(j,:,:) = MISSING
        reflectance_L_SURF(j,:,:) = MISSING
        cycle
      end if
        
      SCAT%pe => pe(:,j)
      SCAT%ze => he(:,j)
      SCAT%te => te(:,j) 
   
      ! Loop over channels
      ! ------------------
      do i = 1, nch
        ! set solar flux
        SCAT%Surface%Base%TSIO%FLUX_FACTOR = flux_factor(i)
 
        ! Make sure kernel weights and parameters are defined
        do n = 1, nkernel
          if (IS_MISSING(kernel_wt(n,i,j),MISSING)) then
            radiance_L_SURF(j,i,:) = MISSING
            reflectance_L_SURF(j,i,:) = MISSING
            cycle
          end if
        end do
        do p = 1, nparam
          if (IS_MISSING(param(p,i,j),MISSING)) then
            radiance_L_SURF(j,i,:) = MISSING
            reflectance_L_SURF(j,i,:) = MISSING
            cycle
          end if
        end do

        if ( verbose > 0 ) then
          print*, 'DO MODIS BRDF'
        end if
        call TWOSTREAM_LANDMODIS(SCAT%Surface,solar_zenith(j,:),&
                               sensor_zenith(j,:),relat_azymuth(j,:),&
                               kernel_wt(1,i,j),kernel_wt(2,i,j),kernel_wt(3,i,j),&
                               reshape(param(:,i,j),(/nparam/)),rc)

        if ( rc /= 0 ) return

        SCAT%wavelength = channels(i)
        SCAT%tau => tau(:,i,j)
        SCAT%ssa => ssa(:,i,j)
        SCAT%fmatrix => pmatrix(:,i,j,:,:)
        SCAT%g => g(:,i,j)
        SCAT%alpha => alpha(:,i,j)
        SCAT%rot => ROT(:,j,i)
        SCAT%depol_ratio => depol_ratio(i)
        SCAT%tauI => tauI(:,i,j)
        SCAT%ssaI => ssaI(:,i,j)
        SCAT%gI => gI(:,i,j)
        SCAT%fmatrixI => pmatrixI(:,i,j,:,:)
        SCAT%tauL => tauL(:,i,j)
        SCAT%ssaL => ssaL(:,i,j)
        SCAT%gL => gL(:,i,j)
        SCAT%fmatrixL => pmatrixL(:,i,j,:,:)
         
        call TWOSTREAM_Run (SCAT, output, ier)

        radiance_L_SURF(j,i,:)    = output%radiance
        reflectance_L_SURF(j,i,:) = output%reflectance

        if ( verbose > 0 ) then
          print *, 'My radiance land modis',radiance_L_SURF(j,i,:), reflectance_L_SURF(j,i,:) 
        end if

        if ( ier /= 0 ) then
          radiance_L_SURF(j,i,:) = MISSING
          reflectance_L_SURF(j,i,:) = MISSING
          cycle
        end if

      end do ! end loop over channels
     
      if ( verbose > 0 ) then
        if ( mod(j-1,1000) == 0 ) then
          print *, '<> TWOSTREAM: ', nint(j*100./nobs), '%'
        end if
      end if

    end do ! Loop over obs

  end subroutine TWOSTREAM_BRDF_LandMODIS


end module TWOSTREAM_BRDF_MODIS
