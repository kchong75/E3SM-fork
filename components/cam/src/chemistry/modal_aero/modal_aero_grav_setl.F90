!===============================================================================
! Calculations related to gravitational settling of aerosol particles. 
!===============================================================================
module modal_aero_grav_setl

  use shr_kind_mod,   only: r8 => shr_kind_r8
  use physconst,      only: gravit
  use perf_mod,       only: t_startf, t_stopf

  implicit none
  private

  public :: interstitial_aero_grav_setl_tend
  public :: modal_aero_gravit_settling_velocity

contains

  !=======================================================================================
  ! Purpose: for interstitial aerosols, calculate the mixing ratio tendencies caused by
  ! gravitational settling.
  !=======================================================================================
  subroutine interstitial_aero_grav_setl_tend( state, pbuf, dt, aerdepdryis_grav, vlc_grv_out, ptend )

    use ppgrid,                  only: pcols, pver
    use constituents,            only: pcnst

    use cam_history,             only: outfld
    use physics_buffer,          only: physics_buffer_desc, pbuf_get_field
    use physics_types,           only: physics_state, physics_ptend, physics_ptend_init

    use physconst,               only: rair
    use modal_aero_data,         only: ntot_amode, nspec_amode
    use constituents,            only: cnst_name
    use modal_aero_data,         only: numptr_amode, lmassptr_amode
    use modal_aero_data,         only: alnsg_amode, sigmag_amode
    use aero_model,              only: drydep_lq, dgnumwet_idx, nmodes, wetdens_ap_idx

    use modal_aero_drydep_utils, only: sedimentation_solver_for_1_tracer, aer_rmax

    ! Arguments

    type(physics_state),target,intent(in) :: state     ! contains atmospheric conditions and aerosol mixing ratios 
    type(physics_buffer_desc), pointer    :: pbuf(:)   ! contains wet diameter and density of interstitial aerosol particles
    real(r8),               intent(in)    :: dt        ! timestep [s]

    real(r8),               intent(out)   :: aerdepdryis_grav(pcols,pcnst)  ! surface deposition flux of 
                                                                            ! interstitial aerosols, [kg/m2/s] or [1/m2/s]

    type(physics_ptend),    intent(out)   :: ptend     ! contains mixing ratio tendencies of interstitial aerosols

    ! Deposition velocities. The last dimension (size = 4) corresponds to the
    ! two attachment states and two moments:
    !   1 - interstitial aerosol, 0th moment (i.e., number)
    !   2 - interstitial aerosol, 3rd moment (i.e., volume/mass)
    !   3 - cloud-borne aerosol,  0th moment (i.e., number)
    !   4 - cloud-borne aerosol,  3rd moment (i.e., volume/mass)
    !
    ! Argument vlc_grv_out has intent(inout) because this subroutine only
    ! calculates values for interstitial aerosols, and we do not want to
    ! touch the velocities of cloud-borne aerosols.

    real(r8),intent(inout) :: vlc_grv_out(pcols,pver,4,ntot_amode)     ! dep velocity of gravitational settling [m/s]

    ! Local variables

    real(r8) :: vlc_grv(pcols,pver,4)     ! dep velocity of gravitational settling [m/s], for  a single mode and moment

    real(r8), pointer :: tair(:,:)   ! air temperture [k]
    real(r8), pointer :: pmid(:,:)   ! air pressure at layer midpoint [Pa]
    real(r8), pointer :: pint(:,:)   ! air pressure at layer interface [Pa]
    real(r8), pointer :: pdel(:,:)   ! layer thickness [Pa]

    integer :: lchnk   ! chunk identifier
    integer :: ncol    ! number of active atmospheric columns
    integer :: lspec   ! index for aerosol number / chem-mass / water-mass
    integer :: imode   ! aerosol mode index
    integer :: icnst   ! tracer index
    integer :: imnt    ! moment of the aerosol size distribution. 0 = number; 3 = volume
    integer :: jvlc    ! index for last dimension of vlc_xxx arrays

    real(r8) ::      rho(pcols,pver) ! air density [kg/m3]
    real(r8) ::     sflx(pcols)      ! surface deposition flux of a single species [kg/m2/s] or [1/m2/s]
    real(r8) :: dqdt_tmp(pcols,pver) ! temporary array to hold tendency for 1 species, [kg/kg/s] or [1/kg/s]

    real(r8), pointer :: dgncur_awet(:,:,:) ! geometric mean wet diameter for number distribution [m]
    real(r8), pointer :: wetdens(:,:,:)     ! wet density of interstitial aerosol [kg/m3]

    real(r8) ::   rad_aer(pcols,pver)  ! volume mean wet radius of interstitial aerosols [m]
    real(r8) ::  dens_aer(pcols,pver)  ! wet density of interstitial aerosols [kg/m3]
    real(r8) ::    sg_aer(pcols,pver)  ! assumed geometric standard deviation of particle size distribution

    real(r8), pointer :: qq(:,:)       ! mixing ratio of a single tracer [kg/kg] or [1/kg]

    character(len=3) :: str

    !-----------------------------------------
    ! Get input variables from state and pbuf
    !-----------------------------------------
    lchnk = state%lchnk
    ncol  = state%ncol

    tair => state%t(:,:)
    pmid => state%pmid(:,:)
    pint => state%pint(:,:)
    pdel => state%pdel(:,:)

    rho(:ncol,:)=  pmid(:ncol,:)/(rair*tair(:ncol,:))

    call pbuf_get_field(pbuf, dgnumwet_idx,   dgncur_awet, start=(/1,1,1/), kount=(/pcols,pver,nmodes/) )
    call pbuf_get_field(pbuf, wetdens_ap_idx, wetdens,     start=(/1,1,1/), kount=(/pcols,pver,nmodes/) )

    !-------------------
    ! Initialize ptend
    !-------------------
    call physics_ptend_init(ptend, state%psetcols, 'interstitial_aero_grav_setl', lq=drydep_lq)

    !-------------------------------------
    ! Loop over aerosol modes
    !-------------------------------------
    do imode = 1, ntot_amode

       !-------------------------------------------------------------------------------------------------
       ! Calculate gravitational settling velocities for interstitial aerosol particles in this mode. 
       ! Note: in each grid cell for each lognormal mode, we have
       !  - one velocity for the number mixing ratio of the mode;
       !  - one velocity for ALL mass mixing ratios in the mode.
       !-----------------------------------------------------------------
       call t_startf('vlc_grv_interstitial')

        rad_aer(:ncol,:) = 0.5_r8*dgncur_awet(:ncol,:,imode) *exp(1.5_r8*(alnsg_amode(imode)**2))
       dens_aer(:ncol,:) = wetdens(:ncol,:,imode)
         sg_aer(:ncol,:) = sigmag_amode(imode)

       ! settling velocity for number mixing ratio

       jvlc = 1  ; imnt = 0
       call modal_aero_gravit_settling_velocity( imnt, ncol, pcols, pver, aer_rmax,     &! in
                                                 tair, pmid, rad_aer, dens_aer, sg_aer, &! in
                                                 vlc_grv(:,:,jvlc)                      )! out

       vlc_grv_out(:ncol,:,jvlc,imode) = vlc_grv(:ncol,:,jvlc)

       ! settling velocity for mass mixing ratios

       jvlc = 2  ; imnt = 3
       call modal_aero_gravit_settling_velocity( imnt, ncol, pcols, pver, aer_rmax,     &! in
                                                 tair, pmid, rad_aer, dens_aer, sg_aer, &! in
                                                 vlc_grv(:,:,jvlc)                      )! out

       vlc_grv_out(:ncol,:,jvlc,imode) = vlc_grv(:ncol,:,jvlc)

       ! Write out those settling velocities

       write(str,'(i0)') imode
       call outfld( 'num_a'//trim(adjustl(str))//'_GVV', vlc_grv(:,:,1), pcols, lchnk )
       call outfld( 'mss_a'//trim(adjustl(str))//'_GVV', vlc_grv(:,:,2), pcols, lchnk )

       call t_stopf('vlc_grv_interstitial')

       !-----------------------------------------------------------
       ! Loop over number + mass species in the mode.
       ! Calculate drydep-induced tendencies; save to ptend.
       !-----------------------------------------------------------
       call t_startf('grav_setl_solve_interstitial')

       do lspec = 0, nspec_amode(imode)

          if (lspec == 0) then  ! number
             icnst = numptr_amode(imode) ; jvlc = 1
          else ! aerosol mass
             icnst = lmassptr_amode(lspec,imode) ; jvlc = 2
          endif

          qq => state%q(:,:,icnst)
          call sedimentation_solver_for_1_tracer( ncol, dt, vlc_grv(:,:,jvlc), qq,    &! in
                                                  rho, tair, pint, pmid, pdel,        &! in
                                                  dqdt_tmp, sflx                      )! out

          aerdepdryis_grav(:ncol,icnst) = sflx(:ncol)
          ptend%lq(icnst) = .true.
          ptend%q(:ncol,:,icnst) = dqdt_tmp(:ncol,:)

          call outfld( trim(cnst_name(icnst))//'GVF',    sflx(:),      pcols, lchnk )
          call outfld( trim(cnst_name(icnst))//'DTQ_GV', dqdt_tmp(:,:),pcols, lchnk )

       enddo ! lspec = 1, nspec_amode(m)

       call t_stopf('grav_setl_solve_interstitial')

    enddo    ! imode = 1, ntot_amode

  end subroutine interstitial_aero_grav_setl_tend

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
