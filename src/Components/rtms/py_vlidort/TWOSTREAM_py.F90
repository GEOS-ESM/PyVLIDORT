!
!  Simple fortran wrapper for the Python interface to 2STREAM 
!  
!
!  Patricia Castellanos
!  April 2025
!.............................................................................


subroutine TWOSTREAM_BRDF_RTLS(km, nch, nobs, channels, plane_parallel, nkernel,nparam, &
                     ROT, depol, alpha, tau, ssa, g, pe, he, te, kernel_wt, param, &
                     solar_zenith, relat_azymuth, sensor_zenith, flux_factor, &
                     MISSING,verbose, radiance_TS_SURF,reflectance_TS_SURF, rc)

    use TWOSTREAM_BRDF_MODIS, only: TWOSTREAM_BRDF_LandMODIS  
    implicit None

  ! !INPUT PARAMETERS:

    integer,          intent(in)            :: km    ! number of vertical levels 
    integer,          intent(in)            :: nch   ! number of channels
    integer,          intent(in)            :: nobs  ! number of observations

    logical,          intent(in)            :: plane_parallel ! do plane parallel flag

    integer,          intent(in)            :: nkernel ! number of kernels
    integer,          intent(in)            :: nparam  ! number of kernel parameters
                    
    real*8,           intent(in)            :: channels(nch)    ! wavelengths [nm]

  !                                                   ! --- Rayleigh Parameters ---
    real*8,           intent(in)            :: ROT(km,nobs,nch) ! rayleigh optical thickness
    real*8,           intent(in)            :: depol(nch)       ! rayleigh depolarization ratio
    real*8,           intent(in)            :: alpha(km,nch,nobs) ! trace gas absorption

  !                                                   ! --- Aerosol Optical Properties ---
    real*8,           intent(in)            :: tau(km,nch,nobs) ! aerosol optical depth
    real*8,           intent(in)            :: ssa(km,nch,nobs) ! single scattering albedo    
    real*8,           intent(in)            :: g(km,nch,nobs)   ! assymetry parameter

    real*8,           intent(in)            :: MISSING          ! MISSING VALUE
    real*8,           intent(in)            :: pe(km+1,nobs)    ! pressure at layer edges [Pa]
    real*8,           intent(in)            :: he(km+1,nobs)    ! height above sea-level  [m]
    real*8,           intent(in)            :: te(km+1,nobs)    ! temperature at layer edges [K]

    real*8,           intent(in)            :: kernel_wt(nkernel,nch,nobs)    ! kernel weights (/fiso,fgeo,fvol/)
    real*8,           intent(in)            :: param(nparam,nch,nobs)         ! Li-Sparse parameters 
                                                                              ! param1 = crown relative height (h/b)
                                                                              ! param2 = shape parameter (b/r)
                         
    real*8,           intent(in)            :: solar_zenith(nobs)  
    real*8,           intent(in)            :: relat_azymuth(nobs) 
    real*8,           intent(in)            :: sensor_zenith(nobs) 

    real*8,           intent(in)            :: flux_factor(nch,nobs) ! solar flux (F0)

    integer,          intent(in)            :: verbose

  ! !OUTPUT PARAMETERS:
    real*8,           intent(out)           :: radiance_TS_SURF(nobs,nch)     ! TOA normalized radiance from VLIDORT using surface module
    real*8,           intent(out)           :: reflectance_TS_SURF(nobs, nch) ! TOA reflectance from VLIDORT using surface module
    integer,          intent(out)           :: rc                             ! return code

    call TWOSTREAM_BRDF_LandMODIS (km, nch, nobs, channels, plane_parallel, &
                                   ROT, depol, alpha, tau, ssa, g, pe, he, te, &
                                   kernel_wt, param, &
                                   solar_zenith, &
                                   relat_azymuth, &
                                   sensor_zenith, &
                                   flux_factor, &
                                   MISSING,verbose, &
                                   radiance_TS_SURF, &
                                   reflectance_TS_SURF, &
                                   rc )  


end subroutine TWOSTREAM_BRDF_RTLS

subroutine TWOSTREAM_LAMBERT_DRIVER(km, nch, nobs, channels, plane_parallel, &
                     ROT, depol, alpha, tau, ssa, g, pe, he, te, albedo, &
                     solar_zenith, relat_azymuth, sensor_zenith, flux_factor, &
                     MISSING,verbose, radiance_TS_SURF,reflectance_TS_SURF, rc)

    use TWOSTREAM_LAMBERT, only: TWOSTREAM_Lambert_Surface
    implicit None

  ! !INPUT PARAMETERS:

    integer,          intent(in)            :: km    ! number of vertical levels 
    integer,          intent(in)            :: nch   ! number of channels
    integer,          intent(in)            :: nobs  ! number of observations

    logical,          intent(in)            :: plane_parallel ! do plane parallel flag

    real*8,           intent(in)            :: channels(nch)    ! wavelengths [nm]

!                                                   ! --- Rayleigh Parameters ---
    real*8,           intent(in)            :: ROT(km,nobs,nch) ! rayleigh optical thickness
    real*8,           intent(in)            :: depol(nch)       ! rayleigh depolarization ratio

    real*8,           intent(in)            :: alpha(km,nobs,nch)       ! trace gas absorption
  !                                                   ! --- Aerosol Optical Properties ---
    real*8,           intent(in)            :: tau(km,nch,nobs) ! aerosol optical depth
    real*8,           intent(in)            :: ssa(km,nch,nobs) ! single scattering albedo    
    real*8,           intent(in)            :: g(km,nch,nobs)   ! assymetry parameter

    real*8,           intent(in)            :: MISSING          ! MISSING VALUE
    real*8,           intent(in)            :: pe(km+1,nobs)    ! pressure at layer edges [Pa]
    real*8,           intent(in)            :: he(km+1,nobs)    ! height above sea-level  [m]
    real*8,           intent(in)            :: te(km+1,nobs)    ! temperature at layer edges [K]

    real*8, target,   intent(in)            :: albedo(nobs,nch)       ! surface albedo
                         
    real*8,           intent(in)            :: solar_zenith(nobs)  
    real*8,           intent(in)            :: relat_azymuth(nobs) 
    real*8,           intent(in)            :: sensor_zenith(nobs) 

    real*8,           intent(in)            :: flux_factor(nch,nobs) ! solar flux (F0)

    integer,          intent(in)            :: verbose

  ! !OUTPUT PARAMETERS:
    real*8,           intent(out)           :: radiance_TS_SURF(nobs,nch)     ! TOA normalized radiance from TWOSTREAM using surface module
    real*8,           intent(out)           :: reflectance_TS_SURF(nobs, nch) ! TOA reflectance from TWOSTREAM using surface module
    integer,          intent(out)           :: rc                             ! return code

    call TWOSTREAM_Lambert_Surface (km, nch, nobs, channels, plane_parallel, &
                                   ROT, depol, alpha, tau, ssa, g, pe, he, te, &
                                   albedo, &
                                   solar_zenith, &
                                   relat_azymuth, &
                                   sensor_zenith, &
                                   flux_factor, &
                                   MISSING,verbose, &
                                   radiance_TS_SURF, &
                                   reflectance_TS_SURF, &
                                   rc )  


end subroutine TWOSTREAM_LAMBERT_DRIVER

