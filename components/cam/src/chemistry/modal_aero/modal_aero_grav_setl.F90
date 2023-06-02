!===============================================================================
! Calculations of 
!===============================================================================
module modal_aero_grav_setl

  use shr_kind_mod,   only: r8 => shr_kind_r8
  use ppgrid,         only: pcols
  use physconst,      only: gravit

  implicit none
  private

  public :: modal_aero_gravit_settling_velocity

contains
  
  !==========================================================================
  ! Calculate particle velocity of gravitational settling
  !==========================================================================
  subroutine modal_aero_gravit_settling_velocity( moment, ncol, pcols, nver, radius_max,           &! in
                                                  tair, pmid, radius_part, density_part, sig_part, &! in
                                                  vlc_grv                                          )! out

    use modal_aero_drydep_utils, only: radius_for_moment, air_dynamic_viscosity, slip_correction_factor

    ! Arguments

    integer,  intent(in) :: moment       ! moment of size distribution (0 for number, 2 for surface area, 3 for volume)
    integer,  intent(in) :: ncol         ! # of grid columns to do calculations for
    integer,  intent(in) :: pcols        ! dimension size (# of columns) 
    integer,  intent(in) :: nver         ! dimension size (# of model layers) 

    real(r8), intent(in) :: radius_max        ! upper bound of radius used for the calculation of deposition velocity [m]

    real(r8), intent(in) :: tair(pcols,nver)    ! air temperature [K]
    real(r8), intent(in) :: pmid(pcols,nver)    ! air pressure [Pa]

    real(r8), intent(in) :: radius_part(pcols,nver)    ! mean (volume or number) particle radius [m]
    real(r8), intent(in) :: density_part(pcols,nver)   ! density of particle material [kg/m3]
    real(r8), intent(in) :: sig_part(pcols,nver)       ! geometric standard deviation of particle size distribution

    real(r8), intent(out) :: vlc_grv(pcols,nver)    ! gravitational deposition velocity [m/s]

    ! Local Variables

    integer  :: ii, kk        ! grid column and layer indices
    real(r8) :: vsc_dyn_atm   ! Dynamic viscosity of air [kg m-1 s-1]
    real(r8) :: slp_crc       ! Slip correction factor [unitless]
    real(r8) :: radius_moment ! median radius for moment [m]

    !-----------
    do kk=1,nver
    do ii=1,ncol

       vsc_dyn_atm = air_dynamic_viscosity( tair(ii,kk) )

       radius_moment = radius_for_moment( moment,sig_part(ii,kk),radius_part(ii,kk),radius_max )

       slp_crc = slip_correction_factor( vsc_dyn_atm, pmid(ii,kk), tair(ii,kk), radius_moment ) 

       vlc_grv(ii,kk) = gravit_settling_velocity( radius_moment, density_part(ii,kk),   &
                                                  slp_crc, vsc_dyn_atm, sig_part(ii,kk) )
    enddo
    enddo
    !----

  end subroutine modal_aero_gravit_settling_velocity

  !=======================================================================================
  ! Calculate the bulk gravitational settling velocity [m s-1] 
  !  - using the terminal velocity of sphere falling in a fluid based on Stokes's law and 
  !  - taking into account the influces of size distribution.
  !=======================================================================================
  real(r8) function gravit_settling_velocity( particle_radius, particle_density,               &
                                              slip_correction, dynamic_viscosity, particle_sig )

    real(r8),intent(in) :: particle_radius   ! [m]
    real(r8),intent(in) :: particle_density  ! [kg/m3]
    real(r8),intent(in) :: slip_correction   ! [unitless]
    real(r8),intent(in) :: dynamic_viscosity ! [kg/m/s]
    real(r8),intent(in) :: particle_sig      ! geometric standard deviation of particle size distribution

    real(r8) :: lnsig         ! ln(particle_sig)
    real(r8) :: dispersion    ! accounts for influence of size dist dispersion on bulk settling velocity
                              ! assuming radius_part is number mode radius * exp(1.5 ln(sigma))

    ! Calculate terminal velocity following, e.g., 
    !  -  Seinfeld and Pandis (1997),  p. 466
    !  - Zhang L. et al. (2001), DOI: 10.1016/S1352-2310(00)00326-5, Eq. 2.

    gravit_settling_velocity = (4.0_r8/18.0_r8) * particle_radius*particle_radius* &
                               particle_density*gravit*slip_correction&
                               /dynamic_viscosity

    ! Account for size distribution (i.e., we are calculating the bulk velocity
    ! for a particle population instead of a single particle).

    lnsig = log(particle_sig)
    dispersion = exp(2._r8*lnsig*lnsig)

    gravit_settling_velocity = gravit_settling_velocity * dispersion

    !---------

  end function gravit_settling_velocity

end module modal_aero_grav_setl
