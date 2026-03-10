module VLIDORT_LAMBERT
!
!  Subroutines to call the VLIDORT lambertian surface preprocessor
!  and then call VLIDORT.s
!
!.............................................................................
  implicit NONE

  PUBLIC VLIDORT_Lambert_pmom
  PUBLIC VLIDORT_Lambert_pmatrix

  contains

  logical function IS_MISSING(x,MISSING)
    real*8, intent(in)     :: x
    real*8, intent(in)     :: MISSING

    IS_MISSING = abs(x/MISSING-1)<0.001
    return
  end function IS_MISSING

  subroutine VLIDORT_Lambert_pmom (km, nch, nobs, ngeom, channels, nstreams, plane_parallel, nMom, &
                     nPol, NSTOKES, ROT, depol, alpha, tau, ssa, pmom, tauI, ssaI, pmomI, tauL, ssaL, pmomL, &
                     pe, he, te, albedo, &
                     solar_zenith, relat_azymuth, sensor_zenith, flux_factor,  &
                     MISSING,verbose,radiance_VL,reflectance_VL, Q, U, rc, &
                     DO_2OS_CORRECTION, DO_BOA)
  !
  ! Uses VLIDORT in scalar mode to compute OMI aerosol TOA radiances.
  !
    use VLIDORT_ScatMod

    implicit NONE

    logical, parameter            :: scalar = .true.   ! is the surface kernel scalar

  ! !INPUT PARAMETERS:

    integer,          intent(in)  :: km    ! number of levels on file
    integer,          intent(in)  :: nch   ! number of channels
    integer,          intent(in)  :: nobs  ! number of observations
    integer,          intent(in)  :: ngeom  ! number of geometries

    logical,          intent(in)  :: plane_parallel ! do plane parallel flag
                                        
    integer, target,  intent(in)  :: nMom             ! number of phase function moments     
    integer, target,  intent(in)  :: nPol  ! number of scattering matrix components                               
    integer,          intent(in)  :: nstreams  ! number of half space streams
    integer,          intent(in)  :: NSTOKES   ! number of stokes vectors. 1 =scalar 3 or 4 =vector

    real*8, target,   intent(in)  :: channels(nch)    ! wavelengths [nm]

!                                                     ! --- Rayleigh Parameters ---
    real*8, target,   intent(in)  :: ROT(km,nobs,nch) ! rayleigh optical thickness
    real*8, target,   intent(in)  :: depol(nch)       ! rayleigh depolarization ratio used in phase matrix

!                                                   ! --- Trace Gas Absorption ---
    real*8, target,   intent(in)  :: alpha(km,nch,nobs) ! trace gas absoprtion optical thickness

  !                                                   ! --- Mie Parameters ---
    real*8, target,   intent(in)  :: tau(km,nch,nobs) ! aerosol optical depth
    real*8, target,   intent(in)  :: ssa(km,nch,nobs) ! single scattering albedo
    real*8, target,   intent(in)  :: pmom(km,nch,nobs,nMom,nPol) !components of the scat phase matrix

    real*8, target,   intent(in)  :: tauI(km,nch,nobs) ! ice cloud optical depth
    real*8, target,   intent(in)  :: ssaI(km,nch,nobs) ! ice cloud single scattering albedo
    real*8, target,   intent(in)  :: pmomI(km,nch,nobs,nMom,nPol) !ice cloud components of the scat phase matrix

    real*8, target,   intent(in)  :: tauL(km,nch,nobs) ! liquid cloud optical depth
    real*8, target,   intent(in)  :: ssaL(km,nch,nobs) ! liquid cloud single scattering albedo
    real*8, target,   intent(in)  :: pmomL(km,nch,nobs,nMom,nPol) !liquid cloud components of the scat phase matrix

    real*8, target,   intent(in)  :: pe(km+1,nobs)    ! pressure at layer edges [Pa]
    real*8, target,   intent(in)  :: he(km+1,nobs)    ! height above sea-level  [m]
    real*8, target,   intent(in)  :: te(km+1,nobs)    ! temperature at layer edges [K]
    
    real*8, target,   intent(in)  :: MISSING          ! MISSING VALUE                                      
    real*8, target,   intent(in)  :: solar_zenith(nobs,ngeom)  
    real*8, target,   intent(in)  :: relat_azymuth(nobs,ngeom) 
    real*8, target,   intent(in)  :: sensor_zenith(nobs,ngeom) 

    real*8, target,   intent(in)  :: albedo(nobs,nch,ngeom)       ! surface albedo

    real*8,           intent(in)  :: flux_factor(nch) ! solar flux (F0)
    
    integer,          intent(in)  :: verbose

  ! !OUTPUT PARAMETERS:

    real*8,           intent(out) :: radiance_VL(nobs,nch,ngeom)       ! TOA normalized radiance from VLIDORT
    integer,          intent(out) :: rc                          ! return code
    real*8,           intent(out) :: reflectance_VL(nobs,nch,ngeom)   ! TOA reflectance from VLIDORT

    real*8,           intent(out) :: Q(nobs, nch, ngeom)                   ! Stokes parameter Q
    real*8,           intent(out) :: U(nobs, nch, ngeom)                   ! Stokes parameter U

  ! !OPTIONAL PARAMETERS
    logical, optional, intent(in) :: DO_2OS_CORRECTION
    logical, optional, intent(in) :: DO_BOA

  !                               ---  
    integer                       :: i,j,ier
   
    type(VLIDORT_scat)            :: SCAT
    type(VLIDORT_output)          :: output  

    rc = 0
    ier = 0

    if (present(DO_BOA)) SCAT%DO_BOA = DO_BOA
    SCAT%Surface%Base%NSTREAMS = nstreams
    SCAT%Surface%Base%DO_PLANE_PARALLEL = plane_parallel
    SCAT%Surface%Base%N_USER_OBSGEOMS   = ngeom
    SCAT%Surface%Base%NBEAMS            = ngeom
    SCAT%Surface%Base%N_USER_STREAMS    = ngeom
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

       ! Make sure albedo and angles are available
       ! -----------------------------------------
       if ( IS_MISSING(solar_zenith(j,1),MISSING)  .OR. & 
            IS_MISSING(sensor_zenith(j,1),MISSING) .OR. &
            IS_MISSING(relat_azymuth(j,1),MISSING)  )  then

          radiance_VL(j,:,:) = MISSING
          reflectance_VL(j,:,:) = MISSING
          Q(j,:,:) = MISSING
          U(j,:,:) = MISSING
          cycle

        end if
        
        SCAT%pe => pe(:,j)
        SCAT%ze => he(:,j)
        SCAT%te => te(:,j) 
   
        ! Loop over channels
        ! ------------------
        do i = 1, nch 
          ! set solar flux
          SCAT%Surface%Base%VIO%VLIDORT_FixIn%SunRays%TS_FLUX_FACTOR = flux_factor(i)         

          ! Mare sure albedo is defined
          ! ---------------------------
          if ( IS_MISSING(albedo(j,i,1),MISSING) ) then
            radiance_VL(j,i,:) = MISSING
            reflectance_VL(j,i,:) = MISSING
            Q(j,:,:) = MISSING
            U(j,:,:) = MISSING
            cycle
          end if

          call VLIDORT_SurfaceLamb(SCAT%Surface,albedo(j,i,:),solar_zenith (j,:),sensor_zenith(j,:),&
                                 relat_azymuth(j,:),scalar)

          SCAT%wavelength = channels(i)
          SCAT%rot => ROT(:,j,i)
          SCAT%depol_ratio => depol(i)
          SCAT%tau => tau(:,i,j)
          SCAT%ssa => ssa(:,i,j)
          SCAT%pmom => pmom(:,i,j,:,:)
          SCAT%tauI => tauI(:,i,j)
          SCAT%ssaI => ssaI(:,i,j)
          SCAT%pmomI => pmomI(:,i,j,:,:)
          SCAT%tauL => tauL(:,i,j)
          SCAT%ssaL => ssaL(:,i,j)
          SCAT%pmomL => pmomL(:,i,j,:,:)

          call VLIDORT_Run (SCAT, output, ier)

          if (SCAT%DO_BOA) then
                radiance_VL(j,i,:)    = output%BOA_radiance
                reflectance_VL(j,i,:) = output%BOA_reflectance
                if (NSTOKES .eq. 1) then
                    ! this should be done upstream,  but impose here explicitly just in case
                    Q(j,i,:)                   = 0.0
                    U(j,i,:)                   = 0.0
                else
                    Q(j,i,:)                   = output%BOA_Q
                    U(j,i,:)                   = output%BOA_U
                end if
          else
                radiance_VL(j,i,:)    = output%radiance
                reflectance_VL(j,i,:) = output%reflectance
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
            radiance_VL(j,i,:) = MISSING
            reflectance_VL(j,i,:) = MISSING
            Q(j,i,:) = MISSING
            U(j,i,:) = MISSING
            cycle
          end if

        end do ! end loop over channels
       
        if ( verbose > 0 ) then
          if ( mod(j-1,1000) == 0 ) then
            print *, '<> VLIDORT Scalar: ', nint(j*100./nobs), '%'
          end if
        end if

    end do ! Loop over obs

  end subroutine VLIDORT_Lambert_pmom

  !..........................................................................

  subroutine VLIDORT_Lambert_pmatrix (km, nch, nobs, ngeom, channels, nstreams, plane_parallel, nAng,  &
                     nPol, NSTOKES, InAngles, ROT, depol, alpha, tau, ssa, pmatrix, tauI, ssaI, pmatrixI, tauL, ssaL, pmatrixL, &
                     pe, he, te, albedo, &
                     solar_zenith, relat_azymuth, sensor_zenith, flux_factor, &
                     MISSING,verbose, radiance_VL, reflectance_VL, Q ,U, rc, &
                     DO_2OS_CORRECTION, DO_BOA)
  !
  ! Place holder.
  !
     use VLIDORT_ScatMod
   
     implicit NONE

    logical, parameter            :: scalar = .true.      ! is the surface BRDF kernel scalar

  ! !INPUT PARAMETERS:

    integer,          intent(in)  :: km    ! number of levels on file
    integer,          intent(in)  :: nch   ! number of channels
    integer,          intent(in)  :: nobs  ! number of observations
    integer,          intent(in)  :: ngeom  ! number of geometries

    logical,          intent(in)  :: plane_parallel ! do plane parallel flag

    integer, target,  intent(in)  :: nAng  ! number of phase function angles
    integer, target,  intent(in)  :: nPol  ! number of components                               
    integer,          intent(in)  :: nstreams  ! number of half space streams
    integer,          intent(in)  :: NSTOKES   ! number of stokes vectors. 1 =scalar 3 or 4 =vector
                    
    real*8, target,   intent(in)  :: channels(nch)    ! wavelengths [nm]
    real*8, target,   intent(in)  :: InAngles(nAng)   ! phase matrix angles

!                                                     ! --- Rayleigh Parameters ---
    real*8, target,   intent(in)  :: ROT(km,nobs,nch) ! rayleigh optical thickness
    real*8, target,   intent(in)  :: depol(nch)       ! rayleigh depolarization ratio used in phase matrix

!                                                   ! --- Trace Gas Absorption ---
    real*8, target,   intent(in)  :: alpha(km,nch,nobs) ! trace gas absorption

  !                                                   ! --- Mie Parameters ---
    real*8, target,   intent(in)  :: tau(km,nch,nobs) ! aerosol optical depth
    real*8, target,   intent(in)  :: ssa(km,nch,nobs) ! single scattering albedo
    real*8, target,   intent(in)  :: pmatrix(km,nch,nobs,nAng,nPol) !components of the scat phase matrix    

    real*8, target,   intent(in)  :: tauI(km,nch,nobs) ! ice cloud optical depth
    real*8, target,   intent(in)  :: ssaI(km,nch,nobs) ! ice cloud single scattering albedo
    real*8, target,   intent(in)  :: pmatrixI(km,nch,nobs,nAng,nPol) !ice cloud components of the scat phase matrix

    real*8, target,   intent(in)  :: tauL(km,nch,nobs) ! liquid cloud optical depth
    real*8, target,   intent(in)  :: ssaL(km,nch,nobs) ! liquid cloud single scattering albedo
    real*8, target,   intent(in)  :: pmatrixL(km,nch,nobs,nAng,nPol) !liquid cloud components of the scat phase matrix


    real*8, target,   intent(in)  :: MISSING          ! MISSING VALUE
    real*8, target,   intent(in)  :: pe(km+1,nobs)    ! pressure at layer edges [Pa]
    real*8, target,   intent(in)  :: he(km+1,nobs)    ! height above sea-level  [m]
    real*8, target,   intent(in)  :: te(km+1,nobs)    ! temperature at layer edges [K]

                         
    real*8, target,   intent(in)  :: solar_zenith(nobs,ngeom)  
    real*8, target,   intent(in)  :: relat_azymuth(nobs,ngeom) 
    real*8, target,   intent(in)  :: sensor_zenith(nobs,ngeom) 

    real*8,           intent(in)  :: flux_factor(nch) ! solar flux (F0)
    real*8, target,   intent(in)  :: albedo(nobs,nch,ngeom)       ! surface albedo

    integer,          intent(in)  :: verbose

  ! !OUTPUT PARAMETERS:

    real*8,           intent(out) :: radiance_VL(nobs,nch,ngeom)       ! TOA radiance from VLIDORT
    integer,          intent(out) :: rc                          ! return code
    real*8,           intent(out) :: reflectance_VL(nobs, nch, ngeom)   ! TOA reflectance from VLIDORT
    real*8,           intent(out) :: Q(nobs, nch, ngeom)   ! Q Stokes component
    real*8,           intent(out) :: U(nobs, nch, ngeom)   ! U Stokes component

  ! !OPTIONAL PARAMETERS
    logical, optional, intent(in) :: DO_2OS_CORRECTION 
    logical, optional, intent(in) :: DO_BOA
    
  !                               ---
    
    integer             :: i,j, ier 
    
    type(VLIDORT_scat) :: SCAT
    type(VLIDORT_output)  :: output  

  
    rc = 0
    ier = 0
    if (present(DO_BOA)) SCAT%DO_BOA = DO_BOA
    SCAT%Surface%Base%NSTREAMS = nstreams
    SCAT%Surface%Base%DO_PLANE_PARALLEL = plane_parallel
    SCAT%Surface%Base%N_USER_OBSGEOMS   = ngeom
    SCAT%Surface%Base%NBEAMS            = ngeom
    SCAT%Surface%Base%N_USER_STREAMS   = ngeom
    SCAT%Surface%Base%N_USER_RELAZMS    = ngeom
    SCAT%Surface%Base%Max_InAngles      = nAng
    SCAT%Surface%Base%USEFMAT = .true.
    call VLIDORT_Init( SCAT%Surface%Base, km, rc, SCAT%DO_BOA)
    if ( rc /= 0 ) return

    SCAT%N_InAngles    = nAng
    SCAT%InAngles => InAngles
    SCAT%nMom    = SCAT%Surface%Base%NGREEK_MOMENTS_INPUT  ! using the default value
    SCAT%nPol = nPol
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

    do j = 1, nobs
       
       ! Make sure albedo and angles are available
       ! -----------------------------------------
       if ( IS_MISSING(solar_zenith(j,1),MISSING)  .OR. & 
            IS_MISSING(sensor_zenith(j,1),MISSING) .OR. &
            IS_MISSING(relat_azymuth(j,1),MISSING)  )  then

          radiance_VL(j,:,:) = MISSING
          reflectance_VL(j,:,:) = MISSING
          Q(j,:,:) = MISSING
          U(j,:,:) = MISSING
          cycle

        end if
       
       SCAT%pe => pe(:,j)
       SCAT%ze => he(:,j)
       SCAT%te => te(:,j) 

       do i = 1, nch
         ! set solar flux
         SCAT%Surface%Base%VIO%VLIDORT_FixIn%SunRays%TS_FLUX_FACTOR = flux_factor(i)

         if ( IS_MISSING(albedo(j,i,1),MISSING) ) then
                radiance_VL(j,i,:) = MISSING
                reflectance_VL(j,i,:) = MISSING
                Q(j,i,:) = MISSING
                U(j,i,:) = MISSING
                cycle
         end if

         call VLIDORT_SurfaceLamb(SCAT%Surface,albedo(j,i,:),solar_zenith (j,:),sensor_zenith(j,:),&
                                 relat_azymuth(j,:),scalar)

          SCAT%wavelength = channels(i)  
          SCAT%rot => ROT(:,j,i)   
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

          call VLIDORT_Run (SCAT, output, ier)

          if (SCAT%DO_BOA) then
            radiance_VL(j,i,:)         = output%BOA_radiance
            reflectance_VL(j,i,:)      = output%BOA_reflectance
            Q(j,i,:)                   = output%BOA_Q
            U(j,i,:)                   = output%BOA_U                          
          else
            radiance_VL(j,i,:)         = output%radiance
            reflectance_VL(j,i,:)      = output%reflectance
            Q(j,i,:)                   = output%Q
            U(j,i,:)                   = output%U                
          end if

          if ( ier /= 0 ) then
            radiance_VL(j,i,:) = MISSING
            reflectance_VL(j,i,:) = MISSING
            Q(j,i,:) = MISSING
            U(j,i,:) = MISSING
            cycle
          end if

          end do ! end loop over channels
       
          if ( verbose > 0 ) then
             if ( mod(j-1,1000) == 0 ) then
                print *, '<> VLIDORT Vector: ', nint(j*100./nobs), '%'
             end if
          end if

    end do ! Loop over obs

  end subroutine VLIDORT_Lambert_pmatrix

  !.............................................................................


end module VLIDORT_LAMBERT
