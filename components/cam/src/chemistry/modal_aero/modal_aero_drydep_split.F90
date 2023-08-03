!===============================================================================
! Calculations of 
!  - changes in interstitial aerosol mixing ratios caused by
!    gravitational settling and turbulent dry deposition of aerosol particles,
!  - changes in cloud-borne aerosol mixing ratios caused by 
!    gravitational settling and turbulent dry deposition of cloud droplets.
!===============================================================================
module modal_aero_drydep_split

  use shr_kind_mod,   only: r8 => shr_kind_r8
  use constituents,   only: pcnst
  use ppgrid,         only: pcols, pver, pverp
  use perf_mod,       only: t_startf, t_stopf
  use cam_history,    only: outfld

  implicit none
  private

  public :: modal_aero_drydep_register
  public :: aero_model_drydep_interstitial
  public :: aero_model_drydep_cloudborne
  public :: interstitial_aero_grav_setl_tend
  public :: interstitial_aero_turb_dep_velocity


  real(r8),parameter :: radius_max = 50.0e-6_r8

contains

  subroutine modal_aero_drydep_register
  !----------------------------------------------------------------------- 
  ! Purpose: Register pbuf fields for aerosols needed by dry deposition
  ! Author: Hui Wan, July 2023
  !-----------------------------------------------------------------------
  use physics_buffer, only: pbuf_add_field, dtype_r8
  use modal_aero_data,only: ntot_amode

  integer :: idx

  call pbuf_add_field('AMODEGVV',   'global',dtype_r8,(/pcols,pver,4,ntot_amode/), idx)
  call pbuf_add_field('AMODETBV',   'global',dtype_r8,(/pcols,     4,ntot_amode/), idx)
  call pbuf_add_field('AERDEPDRYIS','global',dtype_r8,(/pcols,pcnst/),             idx)
  call pbuf_add_field('AERDEPDRYCW','global',dtype_r8,(/pcols,pcnst/),             idx)

  end subroutine modal_aero_drydep_register

  !=============================================================================
  ! Main subroutine for dry deposition of cloud-borne aerosols. 
  ! Also serves as the interface routine called by EAM's physics driver.
  !=============================================================================
  subroutine aero_model_drydep_cloudborne  ( state, pbuf, ram1, fricvel, dt, aerdepdrycw )

    use physics_types,     only: physics_state
    use physics_buffer,    only: physics_buffer_desc

    use physconst,         only: rair, rhoh2o
    use modal_aero_data,   only: ntot_amode, nspec_amode
    use modal_aero_data,   only: qqcw_get_field, cnst_name_cw
    use modal_aero_data,   only: numptrcw_amode, lmassptrcw_amode

    use modal_aero_drydep_utils, only: sedimentation_solver_for_1_tracer
    use modal_aero_drydep,       only: drydep_diags_for_1_tracer, modal_aero_depvel_part

    ! Arguments

    type(physics_state),target,intent(in) :: state     ! Physics state variables
    type(physics_buffer_desc), pointer    :: pbuf(:)
    real(r8),               intent(in)    :: ram1(pcols)    ! aerodynamical resistance, for turbulent dry deposition velocity [s/m]
    real(r8),               intent(in)    :: fricvel(pcols) ! friction velocity, for  turbulent dry deposition velocity [m/s]
    real(r8),               intent(in)    :: dt        ! time step [s]

    real(r8),intent(out) :: aerdepdrycw(pcols,pcnst)  ! surface deposition flux of cloud-borne  aerosols, [kg/m2/s] or [1/m2/s]

    ! Local variables

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

    real(r8) :: rho(pcols,pver)      ! air density [kg/m3]
    real(r8) :: sflx(pcols)          ! surface deposition flux of a single species [kg/m2/s] or [1/m2/s]
    real(r8) :: dqdt_tmp(pcols,pver) ! temporary array to hold tendency for 1 species, [kg/kg/s] or [1/kg/s]

    real(r8) ::  rad_drop(pcols,pver)  ! cloud droplet radius [m]
    real(r8) :: dens_drop(pcols,pver)  ! cloud droplet density [kg/m3]
    real(r8) ::   sg_drop(pcols,pver)  ! assumed geometric standard deviation of droplet size distribution

    real(r8), pointer :: qq(:,:)            ! mixing ratio of a single tracer [kg/kg] or [1/kg]

    ! Deposition velocities. The last dimension (size = 4) corresponds to the
    ! two attachment states and two moments:
    !   1 - interstitial aerosol, 0th moment (i.e., number)
    !   2 - interstitial aerosol, 3rd moment (i.e., volume/mass)
    !   3 - cloud-borne aerosol,  0th moment (i.e., number)
    !   4 - cloud-borne aerosol,  3rd moment (i.e., volume/mass)

    real(r8) :: vlc_grv(pcols,pver,4)     ! dep velocity of gravitational settling [m/s]
    real(r8)::  vlc_trb(pcols,4)          ! dep velocity of turbulent dry deposition [m/s]
    real(r8) :: vlc_dry(pcols,pver,4)     ! dep velocity, sum of vlc_grv and vlc_trb [m/s]

    !---------------------------
    ! Retrieve input variables
    !---------------------------
    lchnk = state%lchnk
    ncol  = state%ncol

    tair => state%t(:,:)
    pmid => state%pmid(:,:)
    pint => state%pint(:,:)
    pdel => state%pdel(:,:)

    rho(:ncol,:)=  pmid(:ncol,:)/(rair*tair(:ncol,:))

    !---------------------------------------------------------------------------------------
    ! Calculate gravitational settling and dry deposition velocities for cloud droplets 
    ! (and hence the cloud-borne aerosols therein).
    ! There is one set of velocities for number mixing ratios of all aerosol modes
    ! and one set of velocities for all mass mixing ratios of all modes.
    !---------------------------------------------------------------------------------------
    ! *** mean drop radius should eventually be computed from ndrop and qcldwtr

    call t_startf('drydep_vel_cloudborne')

    rad_drop(:,:) = 5.0e-6_r8
    dens_drop(:,:) = rhoh2o
    sg_drop(:,:) = 1.46_r8

    jvlc = 3 ; imnt = 0  ! cloud-borne aerosol number
    call modal_aero_depvel_part( ncol, lchnk, tair, pmid, ram1, fricvel,  &! in
                                 rad_drop, dens_drop, sg_drop, imnt,      &! in
                                 vlc_dry(:,:,jvlc),                       &! out
                                 vlc_trb(:,  jvlc),                       &! out
                                 vlc_grv(:,:,jvlc)                        )! out

    jvlc = 4 ; imnt = 3  ! cloud-borne aerosol volume/mass
    call modal_aero_depvel_part( ncol, lchnk, tair, pmid, ram1, fricvel, &! in
                                 rad_drop, dens_drop, sg_drop, imnt,     &! in
                                 vlc_dry(:,:,jvlc),                      &! out
                                 vlc_trb(:,  jvlc),                      &! out
                                 vlc_grv(:,:,jvlc)                       )! out

    call t_stopf('drydep_vel_cloudborne')

    !----------------------------------------------------------------------------------
    ! Loop over all modes and all aerosol tracers (number + mass species).
    ! Calculate the drydep-induced tendencies, then update the mixing ratios.
    !----------------------------------------------------------------------------------
    call t_startf('drydep_solve_cloudborne')

    do imode = 1, ntot_amode         ! loop over aerosol modes
    do lspec = 0, nspec_amode(imode) ! loop over number + constituents

       if (lspec == 0) then   ! number
           icnst = numptrcw_amode(imode) ; jvlc = 3
       else ! aerosol mass
           icnst = lmassptrcw_amode(lspec,imode) ; jvlc = 4
       endif

       qq => qqcw_get_field(pbuf,icnst,lchnk)
       call sedimentation_solver_for_1_tracer( ncol, dt, vlc_dry(:,:,jvlc), qq,    &! in
                                               rho, tair, pint, pmid, pdel,        &! in
                                               dqdt_tmp, sflx                      )! out

       aerdepdrycw(:ncol,icnst) = sflx(:ncol)

       ! Update mixing ratios here. Recall that mixing ratios of cloud-borne aerosols
       ! are stored in pbuf, not as part of the state variable

       qq(1:ncol,:) = qq(1:ncol,:) + dqdt_tmp(1:ncol,:) * dt

       ! Get and save diagnostics

       call drydep_diags_for_1_tracer( lchnk, ncol, trim(cnst_name_cw(icnst)), &! in
                                       vlc_dry(:,:,jvlc), vlc_trb(:,jvlc),     &! in
                                       vlc_grv(:,:,jvlc), sflx                 )! in

    enddo ! loop over number + constituents
    enddo ! imode = 1, ntot_amode

    call t_stopf('drydep_solve_cloudborne')

  end subroutine aero_model_drydep_cloudborne

  !------------------------------------------------------------
  subroutine aero_model_drydep_interstitial  ( state, pbuf, ram1, fricvel, dt, aerdepdryis, ptend )

    use physics_buffer,          only: physics_buffer_desc, pbuf_get_field
    use physics_types,           only: physics_state, physics_ptend, physics_ptend_init

    use physconst,               only: rair
    use modal_aero_data,         only: ntot_amode, nspec_amode
    use constituents,            only: cnst_name
    use modal_aero_data,         only: numptr_amode, lmassptr_amode
    use modal_aero_data,         only: alnsg_amode, sigmag_amode
    use aero_model,              only: drydep_lq, dgnumwet_idx, nmodes, wetdens_ap_idx
    use modal_aero_drydep_utils, only: sedimentation_solver_for_1_tracer
    use modal_aero_drydep,       only: drydep_diags_for_1_tracer, modal_aero_depvel_part

    ! Arguments

    type(physics_state),target,intent(in) :: state     ! Physics state variables
    type(physics_buffer_desc), pointer    :: pbuf(:)
    real(r8),               intent(in)    :: ram1(pcols)    ! aerodynamical resistance, for turbulent dry deposition velocity [s/m]
    real(r8),               intent(in)    :: fricvel(pcols) ! friction velocity, for  turbulent dry deposition velocity [m/s]
    real(r8),               intent(in)    :: dt        ! time step [s]

    real(r8),               intent(out)   :: aerdepdryis(pcols,pcnst)  ! surface deposition flux of interstitial aerosols, [kg/m2/s] or [1/m2/s]
    type(physics_ptend),    intent(out)   :: ptend     ! indivdual parameterization tendencies

    ! Local variables

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

    real(r8) :: rho(pcols,pver)      ! air density [kg/m3]
    real(r8) :: sflx(pcols)          ! surface deposition flux of a single species [kg/m2/s] or [1/m2/s]
    real(r8) :: dqdt_tmp(pcols,pver) ! temporary array to hold tendency for 1 species, [kg/kg/s] or [1/kg/s]

    real(r8), pointer :: dgncur_awet(:,:,:) ! geometric mean wet diameter for number distribution [m]
    real(r8), pointer :: wetdens(:,:,:)     ! wet density of interstitial aerosol [kg/m3]

    real(r8) ::   rad_aer(pcols,pver)  ! volume mean wet radius of interstitial aerosols [m]
    real(r8) ::  dens_aer(pcols,pver)  ! wet density of interstitial aerosols [kg/m3]
    real(r8) ::    sg_aer(pcols,pver)  ! assumed geometric standard deviation of particle size distribution

    real(r8), pointer :: qq(:,:)            ! mixing ratio of a single tracer [kg/kg] or [1/kg]

    ! Deposition velocities. The last dimension (size = 4) corresponds to the
    ! two attachment states and two moments:
    !   1 - interstitial aerosol, 0th moment (i.e., number)
    !   2 - interstitial aerosol, 3rd moment (i.e., volume/mass)
    !   3 - cloud-borne aerosol,  0th moment (i.e., number)
    !   4 - cloud-borne aerosol,  3rd moment (i.e., volume/mass)

    real(r8) :: vlc_grv(pcols,pver,4)     ! dep velocity of gravitational settling [m/s]
    real(r8)::  vlc_trb(pcols,4)          ! dep velocity of turbulent dry deposition [m/s]
    real(r8) :: vlc_dry(pcols,pver,4)     ! dep velocity, sum of vlc_grv and vlc_trb [m/s]

    character(len=3) :: str
    !--------------------------
    ! Retrieve input variables
    !--------------------------
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
    call physics_ptend_init(ptend, state%psetcols, 'aero_model_drydep_ma', lq=drydep_lq)

    !=====================
    ! interstial aerosols
    !=====================
    do imode = 1, ntot_amode   ! loop over aerosol modes

       !-----------------------------------------------------------------
       ! Calculate gravitational settling and dry deposition velocities for 
       ! interstitial aerosol particles in a single lognormal mode. Note:
       !  One set of velocities for number mixing ratio of the mode;
       !  One set of velocities for all mass mixing ratios of the mode.
       !-----------------------------------------------------------------
       call t_startf('drydep_vel_interstitial')

        rad_aer(1:ncol,:) = 0.5_r8*dgncur_awet(1:ncol,:,imode) *exp(1.5_r8*(alnsg_amode(imode)**2))
       dens_aer(1:ncol,:) = wetdens(1:ncol,:,imode)
         sg_aer(1:ncol,:) = sigmag_amode(imode)

       jvlc = 1  ; imnt = 0  ! interstitial aerosol number
       call modal_aero_depvel_part( ncol, lchnk, tair, pmid, ram1, fricvel, &! in
                                    rad_aer, dens_aer, sg_aer, imnt,        &! in
                                    vlc_dry(:,:,jvlc),                      &! out
                                    vlc_trb(:,  jvlc),                      &! out
                                    vlc_grv(:,:,jvlc)                       )! out

       jvlc = 2  ; imnt = 3  ! interstitial aerosol volume/mass
       call modal_aero_depvel_part( ncol, lchnk, tair, pmid, ram1, fricvel, &! in
                                    rad_aer, dens_aer, sg_aer, imnt,        &! in
                                    vlc_dry(:,:,jvlc),                      &! out
                                    vlc_trb(:,  jvlc),                      &! out
                                    vlc_grv(:,:,jvlc)                       )! out

       call t_stopf('drydep_vel_interstitial')
       !-----------------------------------------------------------
       ! Loop over number + mass species of the mode. 
       ! Calculate drydep-induced tendencies; save to ptend.
       !-----------------------------------------------------------
       call t_startf('drydep_solve_interstitial')

       do lspec = 0, nspec_amode(imode)

          if (lspec == 0) then   ! number
             icnst = numptr_amode(imode) ; jvlc = 1
          else ! aerosol mass
             icnst = lmassptr_amode(lspec,imode) ; jvlc = 2
          endif

          qq => state%q(:,:,icnst)
          call sedimentation_solver_for_1_tracer( ncol, dt, vlc_dry(:,:,jvlc), qq,    &! in
                                                  rho, tair, pint, pmid, pdel,        &! in
                                                  dqdt_tmp, sflx                      )! out

          aerdepdryis(:ncol,icnst) = sflx(:ncol)
          ptend%lq(icnst) = .true.
          ptend%q(:ncol,:,icnst) = dqdt_tmp(:ncol,:)

          ! Get and save diagnostics

          call drydep_diags_for_1_tracer( lchnk, ncol, trim(cnst_name(icnst)), &! in
                                          vlc_dry(:,:,jvlc), vlc_trb(:,jvlc),  &! in
                                          vlc_grv(:,:,jvlc), sflx,             &! in
                                          ptend%q(:,:,icnst),pdel(:,:)         )! in

          call outfld( trim(cnst_name(icnst))//'DDV', vlc_dry(:,:,jvlc), pcols, lchnk )

       enddo ! lspec = 1, nspec_amode(m)

       !---------------------------------
       ! Write out deposition velocities
       !---------------------------------
       write(str,'(i0)') imode
       call outfld( 'num_a'//trim(adjustl(str))//'_GVV', vlc_grv(:,:,1), pcols, lchnk )
       call outfld( 'mss_a'//trim(adjustl(str))//'_GVV', vlc_grv(:,:,2), pcols, lchnk )
       call outfld( 'num_a'//trim(adjustl(str))//'_TBV', vlc_trb(:,  1), pcols, lchnk )
       call outfld( 'mss_a'//trim(adjustl(str))//'_TBV', vlc_trb(:,  2), pcols, lchnk )

       call t_stopf('drydep_solve_interstitial')

    enddo    ! imode = 1, ntot_amode

  end subroutine aero_model_drydep_interstitial

  !------------------------------------------------------------
  subroutine interstitial_aero_grav_setl_tend( state, pbuf, dt, aerdepdryis_grav, vlc_grv_out, ptend )

    use cam_history,             only: outfld
    use physics_buffer,          only: physics_buffer_desc, pbuf_get_field
    use physics_types,           only: physics_state, physics_ptend, physics_ptend_init

    use physconst,               only: rair
    use modal_aero_data,         only: ntot_amode, nspec_amode
    use constituents,            only: cnst_name
    use modal_aero_data,         only: numptr_amode, lmassptr_amode
    use modal_aero_data,         only: alnsg_amode, sigmag_amode
    use aero_model,              only: drydep_lq, dgnumwet_idx, nmodes, wetdens_ap_idx

    use modal_aero_grav_setl,    only: modal_aero_gravit_settling_velocity
    use modal_aero_drydep_utils, only: sedimentation_solver_for_1_tracer

    ! Arguments

    type(physics_state),target,intent(in) :: state     ! Physics state variables
    type(physics_buffer_desc), pointer    :: pbuf(:)
    real(r8),               intent(in)    :: dt        ! time step [s]

    real(r8),               intent(out)   :: aerdepdryis_grav(pcols,pcnst)  ! surface deposition flux of interstitial aerosols, [kg/m2/s] or [1/m2/s]
    type(physics_ptend),    intent(out)   :: ptend     ! indivdual parameterization tendencies

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

    real(r8) :: vlc_grv(pcols,pver,4)     ! dep velocity of gravitational settling [m/s]

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

    real(r8) :: rho(pcols,pver)      ! air density [kg/m3]
    real(r8) :: sflx(pcols)          ! surface deposition flux of a single species [kg/m2/s] or [1/m2/s]
    real(r8) :: dqdt_tmp(pcols,pver) ! temporary array to hold tendency for 1 species, [kg/kg/s] or [1/kg/s]

    real(r8), pointer :: dgncur_awet(:,:,:) ! geometric mean wet diameter for number distribution [m]
    real(r8), pointer :: wetdens(:,:,:)     ! wet density of interstitial aerosol [kg/m3]

    real(r8) ::   rad_aer(pcols,pver)  ! volume mean wet radius of interstitial aerosols [m]
    real(r8) ::  dens_aer(pcols,pver)  ! wet density of interstitial aerosols [kg/m3]
    real(r8) ::    sg_aer(pcols,pver)  ! assumed geometric standard deviation of particle size distribution

    real(r8), pointer :: qq(:,:)            ! mixing ratio of a single tracer [kg/kg] or [1/kg]

    character(len=3) :: str
    !--------------------------
    ! Retrieve input variables
    !--------------------------
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

    !=====================
    ! interstial aerosols
    !=====================
    do imode = 1, ntot_amode   ! loop over aerosol modes

       !-----------------------------------------------------------------
       ! Calculate gravitational settling and dry deposition velocities for 
       ! interstitial aerosol particles in a single lognormal mode. Note:
       !  One set of velocities for number mixing ratio of the mode;
       !  One set of velocities for all mass mixing ratios of the mode.
       !-----------------------------------------------------------------
       call t_startf('vlc_grv_interstitial')

        rad_aer(:ncol,:) = 0.5_r8*dgncur_awet(:ncol,:,imode) *exp(1.5_r8*(alnsg_amode(imode)**2))
       dens_aer(:ncol,:) = wetdens(:ncol,:,imode)
         sg_aer(:ncol,:) = sigmag_amode(imode)

       ! interstitial aerosol number

       jvlc = 1  ; imnt = 0
       call modal_aero_gravit_settling_velocity( imnt, ncol, pcols, pver, radius_max,   &! in
                                                 tair, pmid, rad_aer, dens_aer, sg_aer, &! in
                                                 vlc_grv(:,:,jvlc)                      )! out

       vlc_grv_out(:ncol,:,jvlc,imode) = vlc_grv(:ncol,:,jvlc)

       ! interstitial aerosol volume/mass

       jvlc = 2  ; imnt = 3
       call modal_aero_gravit_settling_velocity( imnt, ncol, pcols, pver, radius_max,   &! in
                                                 tair, pmid, rad_aer, dens_aer, sg_aer, &! in
                                                 vlc_grv(:,:,jvlc)                      )! out

       vlc_grv_out(:ncol,:,jvlc,imode) = vlc_grv(:ncol,:,jvlc)

       call t_stopf('vlc_grv_interstitial')
       !-----------------------------------------------------------
       ! Loop over number + mass species of the mode. 
       ! Calculate drydep-induced tendencies; save to ptend.
       !-----------------------------------------------------------
       call t_startf('grav_setl_solve_interstitial')

       do lspec = 0, nspec_amode(imode)

          if (lspec == 0) then   ! number
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

       !-------------------------------
       ! Write out settling velocities
       !-------------------------------
       write(str,'(i0)') imode
       call outfld( 'num_a'//trim(adjustl(str))//'_GVV', vlc_grv(:,:,1), pcols, lchnk )
       call outfld( 'mss_a'//trim(adjustl(str))//'_GVV', vlc_grv(:,:,2), pcols, lchnk )

       call t_stopf('grav_setl_solve_interstitial')

    enddo    ! imode = 1, ntot_amode

  end subroutine interstitial_aero_grav_setl_tend

  !------------------------------------------------
  subroutine interstitial_aero_turb_dep_velocity( state, pbuf, ram1, fricvel, vlc_grv, vlc_trb )

    use physics_types,    only: physics_state
    use physics_buffer,   only: physics_buffer_desc, pbuf_get_field
    use modal_aero_data,  only: ntot_amode
    use modal_aero_data,  only: alnsg_amode, sigmag_amode
    use aero_model,       only: dgnumwet_idx, wetdens_ap_idx, nmodes

    use modal_aero_turb_drydep, only: modal_aero_turb_drydep_velocity

    type(physics_state),    intent(in) :: state
    type(physics_buffer_desc), pointer :: pbuf(:)

    real(r8),intent(in)  ::    ram1(pcols)     ! aerodynamical resistance used in the calculaiton of turbulent dry deposition velocity [s/m]
    real(r8),intent(in)  :: fricvel(pcols)     ! friction velocity used in the calculaiton of turbulent dry deposition velocity [m/s]

    ! Deposition velocities. The last dimension (size = 4) corresponds to the
    ! two attachment states and two moments:
    !   1 - interstitial aerosol, 0th moment (i.e., number)
    !   2 - interstitial aerosol, 3rd moment (i.e., volume/mass)
    !   3 - cloud-borne aerosol,  0th moment (i.e., number)
    !   4 - cloud-borne aerosol,  3rd moment (i.e., volume/mass)
    !
    ! Argument vlc_trb has intent(inout) because this subroutine only
    ! calculates values for interstitial aerosols, and we do not want to
    ! touch the velocities of cloud-borne aerosols.

    real(r8),intent(in)    :: vlc_grv(pcols,4,ntot_amode)  !  gravitational settling   velocity in lowest layer 
    real(r8),intent(inout) :: vlc_trb(pcols,4,ntot_amode)  !  turbulent dry deposition velocity in lowest layer 

    ! Local variables

    real(r8) :: vlc_dry(pcols)  !  total dry deposition velocity in lowest layer 


    real(r8), pointer :: dgncur_awet(:,:,:) ! geometric mean wet diameter for number distribution [m]
    real(r8), pointer :: wetdens(:,:,:)     ! wet density of interstitial aerosol [kg/m3]

    real(r8) ::   rad_aer(pcols)  ! volume mean wet radius of interstitial aerosols [m]
    real(r8) ::  dens_aer(pcols)  ! wet density of interstitial aerosols [kg/m3]
    real(r8) ::    sg_aer(pcols)  ! assumed geometric standard deviation of particle size distribution

    real(r8) ::    tair(pcols)
    real(r8) ::    pmid(pcols)
    real(r8) ::    pdel(pcols)

    integer :: lchnk   ! chunk identifier
    integer :: ncol    ! number of active atmospheric columns
    integer :: imode   ! aerosol mode index
    integer :: imnt    ! moment of the aerosol size distribution. 0 = number; 3 = volume
    integer :: jvlc    ! index for last dimension of vlc_xxx arrays

    character(len=3) :: str

    lchnk = state%lchnk
    ncol  = state%ncol

    tair  = state%t(:,pver)
    pmid  = state%pmid(:,pver)
    pdel  = state%pdel(:,pver)

    call pbuf_get_field(pbuf, dgnumwet_idx,   dgncur_awet, start=(/1,1,1/), kount=(/pcols,pver,nmodes/) ) 
    call pbuf_get_field(pbuf, wetdens_ap_idx, wetdens,     start=(/1,1,1/), kount=(/pcols,pver,nmodes/) ) 

    !------------------------------------------------------------------------
    ! Calculate turbulent dry deposition velocities of interstitial aerosols.
    ! Note that within a grid cell, for each lognormal mode of the 
    ! aerosol size distribution, we have
    !  - one velocity for number mixing ratio of the mode;
    !  - one velocity for all mass mixing ratios of the mode.
    !-----------------------------------------------------------------
    call t_startf('vlc_trb_interstitial')

    do imode = 1, ntot_amode   ! loop over aerosol modes

        rad_aer(1:ncol) = 0.5_r8*dgncur_awet(1:ncol,pver,imode) *exp(1.5_r8*(alnsg_amode(imode)**2))
       dens_aer(1:ncol) = wetdens(1:ncol,pver,imode)
         sg_aer(1:ncol) = sigmag_amode(imode)

       jvlc = 1  ; imnt = 0  ! interstitial aerosol number
       call modal_aero_turb_drydep_velocity( imnt, ncol, pcols, lchnk, radius_max,  &! in
                                             tair, pmid,                            &! in
                                             rad_aer, dens_aer, sg_aer,             &! in
                                             fricvel, ram1,                         &! in
                                             vlc_grv(:,jvlc,imode), &! in
                                             vlc_trb(:,jvlc,imode), &! out
                                             vlc_dry(:)             )! out

       jvlc = 2  ; imnt = 3  ! interstitial aerosol volume/mass
       call modal_aero_turb_drydep_velocity( imnt, ncol, pcols, lchnk, radius_max,  &! in
                                             tair, pmid,                            &! in
                                             rad_aer, dens_aer, sg_aer,             &! in
                                             fricvel, ram1,                         &! in
                                             vlc_grv(:,jvlc,imode), &! in
                                             vlc_trb(:,jvlc,imode), &! out
                                             vlc_dry(:)             )! out

       !-------------------------------
       ! Write out settling velocities
       !-------------------------------
       write(str,'(i0)') imode
       call outfld( 'num_a'//trim(adjustl(str))//'_TBV', vlc_trb(:,1,imode), pcols, lchnk )
       call outfld( 'mss_a'//trim(adjustl(str))//'_TBV', vlc_trb(:,2,imode), pcols, lchnk )

    enddo ! imode = 1, ntot_amode

    call t_stopf('vlc_trb_interstitial')

  end subroutine interstitial_aero_turb_dep_velocity

end module modal_aero_drydep_split
