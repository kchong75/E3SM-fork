!===============================================================================
!===============================================================================
module modal_aero_drydep_utils

  use shr_kind_mod,   only: r8 => shr_kind_r8
  use ppgrid,         only: pcols, pver, pverp
  use physconst,      only: gravit, rair, pi, boltz
  use perf_mod,       only: t_startf, t_stopf

  implicit none
  public

  real(r8),parameter,public :: aer_rmax = 50.0e-6_r8

contains

  !=======================================================================================
  ! Purpose: Register pbuf fields for aerosols needed by dry deposition
  ! Author: Hui Wan, July 2023
  !=======================================================================================
  subroutine modal_aero_drydep_register

  use physics_buffer, only: pbuf_add_field, dtype_r8
  use modal_aero_data,only: ntot_amode
  use constituents,   only: pcnst

  integer :: idx

  call pbuf_add_field('AMODEGVV',   'global',dtype_r8,(/pcols,pver,4,ntot_amode/), idx)
  call pbuf_add_field('AMODETBV',   'global',dtype_r8,(/pcols,     4,ntot_amode/), idx)
  call pbuf_add_field('AERDEPDRYIS','global',dtype_r8,(/pcols,pcnst/),             idx)
  call pbuf_add_field('AERDEPDRYCW','global',dtype_r8,(/pcols,pcnst/),             idx)

  end subroutine modal_aero_drydep_register


  !=================================================================================
  !=================================================================================
  subroutine outfld_aero_cnst_2d( fld, suffix, lchnk )

  use modal_aero_data,         only: ntot_amode, nspec_amode
  use constituents,            only: cnst_name
  use modal_aero_data,         only: numptr_amode, lmassptr_amode
  use cam_history,             only: outfld
 
  real(r8),        intent(in) :: fld(:,:)
  character(len=*),intent(in) :: suffix
  integer,         intent(in) :: lchnk

  integer :: imode, lspec, icnst

  do imode=1,ntot_amode
  do lspec = 0, nspec_amode(imode)

     if (lspec == 0) then   ! number
        icnst = numptr_amode(imode)
     else ! aerosol mass
        icnst = lmassptr_amode(lspec,imode)
     endif

     call outfld( trim(cnst_name(icnst))//trim(suffix), fld(:,icnst), pcols, lchnk)

  end do
  end do

  end subroutine outfld_aero_cnst_2d  

  !=============================================================================
  ! Numerically solve the sedimentation equation for 1 tracer
  !=============================================================================
  subroutine sedimentation_solver_for_1_tracer( ncol, dt, sed_vel, qq_in,            &! in
                                                rho, tair, pint, pmid, pdel,         &! in
                                                dqdt_sed, sflx                       )! out

    use shr_kind_mod,      only: r8 => shr_kind_r8
    use ppgrid,            only: pcols, pver, pverp
    use mo_spitfire_transport, only: getflx 

    integer , intent(in) :: ncol
    real(r8), intent(in) :: dt
    real(r8), intent(in) :: rho(pcols,pver)           ! air density [kg/m3]
    real(r8), intent(in) :: tair(pcols,pver)          ! air temperature [K]
    real(r8), intent(in) :: pint(pcols,pverp)         ! air pressure at layer interfaces [Pa]
    real(r8), intent(in) :: pmid(pcols,pver)          ! air pressure at layer midpoints  [Pa]
    real(r8), intent(in) :: pdel(pcols,pver)          ! pressure layer thickness [Pa]
    real(r8), intent(in) :: sed_vel(pcols,pver)       ! deposition velocity [m/s]
    real(r8), intent(in) :: qq_in(pcols,pver)         ! tracer mixing ratio, [kg/kg] or [1/kg]

    real(r8), intent(out) :: dqdt_sed(pcols,pver) ! tracer mixing ratio tendency [kg/kg/s] or [1/kg/s]
    real(r8), intent(out) :: sflx(pcols)          ! deposition flux at the Earth's surface [kg/m2/s] or [1/m2/s]

    real (r8), parameter :: mxsedfac = 0.99_r8    ! maximum sedimentation flux factor

    real(r8) :: pvmzaer(pcols,pverp)     ! sedimentation velocity in Pa (positive = down)
    real(r8) :: dtmassflux(pcols,pverp)  ! dt * mass fluxes at layer interfaces (positive = down)

    integer  :: ii,kk

    !---------------------------------------------------------------------------------------
    ! Set sedimentation velocity to zero at the top interface of the model domain.

    pvmzaer(:ncol,1)=0._r8

    ! Assume the sedimentation velocities passed in are velocities
    ! at the bottom interface of each model layer, like an upwind scheme.

    pvmzaer(:ncol,2:pverp) = sed_vel(:ncol,:)

    ! Convert velocity from height coordinate to pressure coordinate; 
    ! units: convert from meters/sec to pascals/sec.
    ! (This was referred to as "Phil's method" in the code before refactoring.)
 
    pvmzaer(:ncol,2:pverp) = pvmzaer(:ncol,2:pverp) * rho(:ncol,:)*gravit

    !------------------------------------------------------
    ! Calculate mass flux * dt at each layer interface
    !------------------------------------------------------
    call getflx(ncol, pint, qq_in, pvmzaer, dt, dtmassflux)

    ! Filter out any negative fluxes from the getflx routine

    do kk = 2,pver
       dtmassflux(:ncol,kk) = max(0._r8, dtmassflux(:ncol,kk))
    enddo

    ! Set values for the upper and lower boundaries

    do ii = 1,ncol
       dtmassflux(ii,1)     = 0                                         ! no flux at model top 
       dtmassflux(ii,pverp) = qq_in(ii,pver) * pvmzaer(ii,pverp) * dt   ! surface flux by upwind scheme
    enddo

    ! Limit the flux out of the bottom of each column:
    ! apply mxsedfac to prevent generating very small negative mixing ratio.
    ! *** Should we include the flux through the top interface, to accommodate thin surface layers?

    do kk = 1,pver
       do ii = 1,ncol
          dtmassflux(ii,kk+1) = min( dtmassflux(ii,kk+1), mxsedfac * qq_in(ii,kk) * pdel(ii,kk))
       enddo
    enddo

    !-----------------------------------------------------------------------
    ! Calculate the mixing ratio tendencies resulting from flux divergence
    !-----------------------------------------------------------------------
    do kk = 1,pver
       do ii = 1,ncol
          dqdt_sed(ii,kk)  = (dtmassflux(ii,kk) - dtmassflux(ii,kk+1)) / (dt * pdel(ii,kk))
       enddo
    enddo

    !-----------------------------------------------------------------------
    ! Convert flux out the bottom to mass units [kg/m2/s]
    !-----------------------------------------------------------------------
    sflx(:ncol) = dtmassflux(:ncol,pverp) / (dt*gravit)

  end subroutine sedimentation_solver_for_1_tracer


  !==============================================================================
  ! Calculate the radius for a moment of a lognormal size distribution
  !==============================================================================
  real(r8) function radius_for_moment( moment,sig_part,radius_part,radius_max)

    integer, intent(in) :: moment        ! moment of the particle size distribution. 0 = number; 3 = volume
    real(r8),intent(in) :: sig_part      ! geometric standard deviation of particle size distribution
    real(r8),intent(in) :: radius_part   ! mean (volume or number) particle radius [m]
    real(r8),intent(in) :: radius_max    ! developer-specified upper bound of mean radius [m]

    real(r8) :: lnsig
    
    lnsig = log(sig_part)
    radius_for_moment = min(radius_max,radius_part)*exp((float(moment)-1.5_r8)*lnsig*lnsig)

  end function radius_for_moment

  !==========================================================================
  ! Calculate dynamic viscosity of air, unit [kg m-1 s-1]. See RoY94 p. 102
  !==========================================================================
  real(r8) function air_dynamic_viscosity( temp )

    real(r8),intent(in) :: temp   ! air temperature [K]

    air_dynamic_viscosity = 1.72e-5_r8 * ( (temp/273.0_r8)**1.5_r8) * 393.0_r8 &
                                          /(temp+120.0_r8)

  end function air_dynamic_viscosity

  !==========================================================================
  ! Calculate kinematic viscosity of air, unit [m2 s-1]
  !==========================================================================
  real(r8) function air_kinematic_viscosity( temp, pres )

    real(r8),intent(in) :: temp   ! air temperature [K]
    real(r8),intent(in) :: pres   ! air pressure [Pa]

    real(r8) :: vsc_dyn_atm  ! dynamic viscosity of air [kg m-1 s-1]
    real(r8) :: rho          ! density of air [kg/m3]

    vsc_dyn_atm = air_dynamic_viscosity( temp)
    rho = pres/rair/temp

    air_kinematic_viscosity = vsc_dyn_atm/rho

  end function air_kinematic_viscosity

  !======================================================
  ! Slip correction factor [unitless]. 
  ! See, e.g., SeP97 p. 464 and Zhang L. et al. (2001),
  ! DOI: 10.1016/S1352-2310(00)00326-5, Eq. (3).
  !======================================================
  real(r8) function slip_correction_factor( dyn_visc, pres, temp, particle_radius ) 

    real(r8),intent(in) :: dyn_visc         ! dynamic viscosity of air [kg m-1 s-1]
    real(r8),intent(in) :: pres             ! air pressure [Pa]
    real(r8),intent(in) :: temp             ! air temperature [K]
    real(r8),intent(in) :: particle_radius  ! particle radius [m]

    real(r8) :: mean_free_path  ! [m]

    mean_free_path = 2.0_r8*dyn_visc/( pres*sqrt( 8.0_r8/(pi*rair*temp) ) )

    slip_correction_factor = 1.0_r8 + mean_free_path * &
                             ( 1.257_r8+0.4_r8*exp(-1.1_r8*particle_radius/mean_free_path) ) / &
                             particle_radius

  end function slip_correction_factor

  !====================================================================
  ! Calculate the Schmidt number of air [unitless], see SeP97 p.972
  !====================================================================
  real(r8) function schmidt_number( temp, pres, radius, vsc_dyn_atm, vsc_knm_atm )

    real(r8),intent(in) :: temp             ! air temperature [K]
    real(r8),intent(in) :: pres             ! air pressure [Pa]
    real(r8),intent(in) :: radius           ! particle radius [m]
    real(r8),intent(in) :: vsc_dyn_atm      ! dynamic viscosity of air [kg m-1 s-1]
    real(r8),intent(in) :: vsc_knm_atm      ! kinematic viscosity of air [m2 s-1]

    real(r8) :: slp_crc   ! slip correction factor [unitless]
    real(r8) :: dff_aer   ! Brownian diffusivity of particle [m2/s], see SeP97 p.474

    slp_crc = slip_correction_factor( vsc_dyn_atm, pres, temp, radius )
    dff_aer = boltz * temp * slp_crc / (6.0_r8*pi*vsc_dyn_atm*radius) 

    schmidt_number = vsc_knm_atm / dff_aer

  end function schmidt_number

end module modal_aero_drydep_utils
