!
!  Simple f77 wrapper for the Python interface to 2STREAM
!
!.............................................................................

subroutine TWOSTREAM_Lambert (km, nch, nobs,channels, nMom, &
                   nPol, tau, ssa, g, pmom, pe, he, te, albedo,            &
                   solar_zenith, relat_azymuth, sensor_zenith, &
                   MISSING,verbose,radiance_L,reflectance_L, ROT, rc)
!
! Uses 2STREAM to compute TOA radiance
!
  use TWOSTREAM_ScatMod

  implicit NONE

! !INPUT PARAMETERS:

  integer,          intent(in)  :: km    ! number of levels on file
  integer,          intent(in)  :: nch   ! number of channels
  integer,          intent(in)  :: nobs  ! number of observations

  real*8, target,   intent(in)  :: channels(nch)    ! wavelengths [nm]

!                                                   ! --- Mie Parameters ---
  real*8, target,   intent(in)  :: tau(km,nch,nobs) ! aerosol optical depth
  real*8, target,   intent(in)  :: ssa(km,nch,nobs) ! single scattering albedo
  real*8, target,   intent(in)  :: g(km,nch,nobs)   ! asymmetry factor

  real*8, target,   intent(in) :: ROT(km,nobs,nch)  ! rayleigh optical thickness
  real*8, target,   intent(in) :: depol_ratio(nch)  ! depolariation ratio

  real*8, target,   intent(in)  :: pe(km+1,nobs)    ! pressure at layer edges [Pa]
  real*8, target,   intent(in)  :: he(km+1,nobs)    ! height above sea-level  [m]
  real*8, target,   intent(in)  :: te(km+1,nobs)    ! temperature at layer edges [K]
  
  real*8, target,   intent(in)  :: MISSING          ! MISSING VALUE                                      
  real*8, target,   intent(in)  :: solar_zenith(nobs)  
  real*8, target,   intent(in)  :: relat_azymuth(nobs) 
  real*8, target,   intent(in)  :: sensor_zenith(nobs) 

  real*8, target,   intent(in)  :: albedo(nobs,nch)       ! surface albedo
  
  integer,          intent(in)  :: verbose


! !OUTPUT PARAMETERS:

  real*8,           intent(out) :: radiance_L(nobs,nch)       ! TOA normalized radiance from LIDORT
  integer,          intent(out) :: rc                          ! return code
  real*8,           intent(out) :: reflectance_L(nobs, nch)   ! TOA reflectance from LIDORT
 
!                               ---  
  integer             :: i,j, ier
 
  type(TWOSTREAM_scat)            :: SCAT
  type(TWOSTREAM_output)          :: output  

#define IS_MISSING(x) (abs(x/MISSING-1)<0.001)
  rc = 0
  ier = 0

  call TWOSTREAM_Init( SCAT%Surface%Base, km, rc)
  if ( rc /= 0 ) then
    write(*,*) 'TWOSTREAM_Init returning with error'
    return
  endif

  do j = 1, nobs

     ! Make sure albedo and angles are available
     ! -----------------------------------------
     if ( IS_MISSING(solar_zenith(j))  .OR. & 
          IS_MISSING(sensor_zenith(j)) .OR. &
          IS_MISSING(relat_azymuth(j))  )  then

        radiance_L(j,:) = MISSING
        reflectance_L(j,:) = MISSING
        cycle

      end if
      
      SCAT%pe => pe(:,j)
      SCAT%ze => he(:,j)
      SCAT%te => te(:,j) 
 
      ! Loop over channels
      ! ------------------
      do i = 1, nch 
       
           ! Mare sure albedo is defined
           ! ---------------------------
           if ( IS_MISSING(albedo(j,i)) ) then
              radiance_L(j,i) = MISSING
              reflectance_L(j,i) = MISSING
              cycle
           end if
           
           call TWOSTREAM_SurfaceLamb(SCAT%Surface,albedo(j,i),solar_zenith (j),sensor_zenith(j),&
                               relat_azymuth(j))
           
           SCAT%wavelength = channels(i)
           SCAT%tau => tau(:,i,j)
           SCAT%ssa => ssa(:,i,j)
           SCAT%g => g(:,i,j)
         
           call TWOSTREAM_Run (SCAT, output, ier)

           radiance_L(j,i)    = output%radiance
           reflectance_L(j,i) = output%reflectance


           if ( ier /= 0 ) then
              radiance_L(j,i) = MISSING
              reflectance_L(j,i) = MISSING
              cycle
           end if

        end do ! end loop over channels
     
        if ( verbose > 0 ) then
           if ( mod(j-1,1000) == 0 ) then
              print *, '<> TWOSTREAM: ', nint(j*100./nobs), '%'
           end if
        end if

  end do ! Loop over obs

end subroutine LIDORT_Scalar_Lambert
!.............................................................................
