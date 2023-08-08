!=========================================================================================
! Calculations of 
!  - changes in interstitial aerosol mixing ratios caused by
!    gravitational settling and turbulent dry deposition of aerosol particles,
!  - changes in cloud-borne aerosol mixing ratios caused by 
!    gravitational settling and turbulent dry deposition of cloud droplets.
!
! This module contains a refactored version of the subroutine aero_model_drydep
! in src/chemistry/modal_aero/aero_model.F90 from 2023, now named aero_model_drydep_main.
! Like in the original code, the calculated dry deposition velocities include effects
! of both turbulent dry deposition and gravitational settling.
! The two processes are numerically solved together using Phil Rasch's SPITFIRE algorithm.
!
! POC of this module: Hui Wan, PNNL.
!=========================================================================================
module modal_aero_drydep_refac

  use shr_kind_mod,   only: r8 => shr_kind_r8
  use constituents,   only: pcnst
  use ppgrid,         only: pcols, pver
  use cam_history,    only: outfld
  use perf_mod,       only: t_startf, t_stopf

  implicit none
  private

  public :: aero_model_drydep_main
  public :: aero_model_drydep_cloudborne
  public :: aero_model_drydep_interstitial

contains

  !=============================================================================
  ! Main subroutine that calculates deposition velocities and solves PDEs
  ! for interstitial and cloud-borne aerosols.
  !=============================================================================
  subroutine aero_model_drydep_main( state, pbuf, cam_in, dt, cam_out, ptend )

    use physics_types,     only: physics_state, physics_ptend
    use physics_buffer,    only: physics_buffer_desc
    use camsrfexch,        only: cam_out_t, cam_in_t

    use aerodep_flx,           only: aerodep_flx_prescribed
    use modal_aero_deposition, only: set_srf_drydep
    use modal_aero_turb_drydep,only: get_gridcell_ram1_fricvel

    type(physics_state),target,intent(in) :: state     ! Physics state variables
    type(physics_buffer_desc), pointer    :: pbuf(:)
    type(cam_in_t),         intent(in)    :: cam_in
    real(r8),               intent(in)    :: dt        ! time step [s]
    type(cam_out_t),        intent(inout) :: cam_out
    type(physics_ptend),    intent(out)   :: ptend     ! indivdual parameterization tendencies

    real(r8) :: fricvel(pcols)     ! friction velocity used in the calculaiton of turbulent dry deposition velocity [m/s]
    real(r8) ::    ram1(pcols)     ! aerodynamical resistance used in the calculaiton of turbulent dry deposition velocity [s/m]

    real(r8) :: aerdepdryis(pcols,pcnst)  ! surface deposition flux of interstitial aerosols, [kg/m2/s] or [1/m2/s]
    real(r8) :: aerdepdrycw(pcols,pcnst)  ! surface deposition flux of cloud-borne  aerosols, [kg/m2/s] or [1/m2/s]

    integer :: lchnk

    lchnk = state%lchnk

    !------------------------------------------------------------------
    ! Calculate/copy ram1 and fricvel for each grid cell in this chunk
    !------------------------------------------------------------------
    call get_gridcell_ram1_fricvel(state, cam_in, ram1, fricvel)
    call outfld( 'RAM1',     ram1(:), pcols, lchnk )
    call outfld( 'airFV', fricvel(:), pcols, lchnk )

    !---------------------------------------------------------------------------------------------
    ! Dry deposition of cloud-borne aerosols.
    ! Note the corresponding mixing ratios are in pbuf, so there is no ptend in the argument list.
    !----------------------------------------------------------------------------------------------
    call aero_model_drydep_cloudborne( state, pbuf, ram1, fricvel, dt, aerdepdrycw )

    !-----------------------------------------------------------------------------------------------------------
    ! Dry deposition of interstitial aerosols.
    ! The corresponding mixing ratios are part of state%q, so ptend is an output variable on the argument list.
    !-----------------------------------------------------------------------------------------------------------
    call aero_model_drydep_interstitial( state, pbuf, ram1, fricvel, dt, aerdepdryis, ptend )

    !------------------------------------------------------------------
    ! Unless the user has specified prescribed aerosol dep fluxes,
    ! copy the fluxes calculated here to cam_out to be passed to other 
    ! components of the Earth System Model.
    !------------------------------------------------------------------
    if (.not.aerodep_flx_prescribed()) then
       call set_srf_drydep( aerdepdryis, aerdepdrycw, cam_out )
    endif

  end subroutine aero_model_drydep_main

  !=============================================================================
  ! Main subroutine for dry deposition of cloud-borne aerosols. 
  !=============================================================================
  subroutine aero_model_drydep_cloudborne  ( state, pbuf, ram1, fricvel, dt, aerdepdrycw )

    use physics_types,     only: physics_state
    use physics_buffer,    only: physics_buffer_desc

    use physconst,         only: rair, rhoh2o
    use modal_aero_data,   only: ntot_amode, nspec_amode
    use modal_aero_data,   only: qqcw_get_field, cnst_name_cw
    use modal_aero_data,   only: numptrcw_amode, lmassptrcw_amode

    use modal_aero_drydep_utils, only: sedimentation_solver_for_1_tracer

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

    real(r8) ::       rho(pcols,pver)  ! air density [kg/m3]
    real(r8) ::      sflx(pcols)       ! surface deposition flux of a single species [kg/m2/s] or [1/m2/s]
    real(r8) ::  dqdt_tmp(pcols,pver)  ! temporary array to hold tendency for 1 species, [kg/kg/s] or [1/kg/s]

    real(r8) ::  rad_drop(pcols,pver)  ! cloud droplet radius [m]
    real(r8) :: dens_drop(pcols,pver)  ! cloud droplet density [kg/m3]
    real(r8) ::   sg_drop(pcols,pver)  ! assumed geometric standard deviation of droplet size distribution

    real(r8), pointer :: qq(:,:)       ! mixing ratio of a single tracer [kg/kg] or [1/kg]

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

       ! Write out deposition velocities

       write(str,'(i0)') imode
       call outfld( 'num_a'//trim(adjustl(str))//'_GVV', vlc_grv(:,:,1), pcols, lchnk )
       call outfld( 'mss_a'//trim(adjustl(str))//'_GVV', vlc_grv(:,:,2), pcols, lchnk )
       call outfld( 'num_a'//trim(adjustl(str))//'_TBV', vlc_trb(:,  1), pcols, lchnk )
       call outfld( 'mss_a'//trim(adjustl(str))//'_TBV', vlc_trb(:,  2), pcols, lchnk )

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

       call t_stopf('drydep_solve_interstitial')

    enddo    ! imode = 1, ntot_amode

  end subroutine aero_model_drydep_interstitial

  !==========================================================================================
  ! Calculate deposition velocities caused by turbulent dry deposition and
  ! gravitational settling of aerosol particles

  ! Reference: 
  !  L. Zhang, S. Gong, J. Padro, and L. Barrie:
  !  A size-seggregated particle dry deposition scheme for an atmospheric aerosol module
  !  Atmospheric Environment, 35, 549-560, 2001.
  !
  ! History: 
  !  - Original version by X. Liu.
  !  - Calculations for gravitational and turbulent dry deposition separated into 
  !    different subroutines by Hui Wan, 2023.
  !==========================================================================================
  subroutine modal_aero_depvel_part( ncol, lchnk, tair, pmid, ram1, fricvel,       &! in
                                     radius_part, density_part, sig_part, moment,  &! in
                                     vlc_dry, vlc_trb, vlc_grv                     )! out

    use modal_aero_grav_setl,   only: modal_aero_gravit_settling_velocity
    use modal_aero_turb_drydep, only: modal_aero_turb_drydep_velocity
    use modal_aero_drydep_utils,only: aer_rmax

    integer,  intent(in) :: moment       ! moment of size distribution (0 for number, 2 for surface area, 3 for volume)
    integer,  intent(in) :: ncol         ! # of grid columns to do calculations for
    integer,  intent(in) :: lchnk        ! chunk index 

    real(r8), intent(in) :: tair(pcols,pver)    ! air temperature [K]
    real(r8), intent(in) :: pmid(pcols,pver)    ! air pressure [Pa]

    real(r8), intent(in) :: ram1(pcols)     ! aerodynamical resistance [s/m]
    real(r8), intent(in) :: fricvel(pcols)  ! sfc friction velocity from land model [m/s]

    real(r8), intent(in) :: radius_part(pcols,pver)    ! mean (volume or number) particle radius [m]
    real(r8), intent(in) :: density_part(pcols,pver)   ! density of particle material [kg/m3]
    real(r8), intent(in) :: sig_part(pcols,pver)       ! geometric standard deviation of particle size distribution

    real(r8), intent(out) :: vlc_grv(pcols,pver)    ! gravitational deposition velocity [m/s]
    real(r8), intent(out) :: vlc_dry(pcols,pver)    ! total dry deposition velocity [m/s]
    real(r8), intent(out) :: vlc_trb(pcols)         ! turbulent dry deposition velocity [m/s]

    ! use a maximum radius of 50 microns when calculating deposition velocity


    !------------------------------------------------------------------------------------
    ! Calculate deposition velocity of gravitational settling in all grid layers
    !------------------------------------------------------------------------------------
    call modal_aero_gravit_settling_velocity( moment, ncol, pcols, pver, aer_rmax,             &! in
                                              tair, pmid, radius_part, density_part, sig_part, &! in
                                              vlc_grv                                          )! out

    ! vlc_dry in layers above the lowest is just the gravitational settling velocity
    vlc_dry(:ncol,1:pver-1) = vlc_grv(:ncol,1:pver-1)

    !------------------------------------------------------------------------------------
    ! For the lowest model layer:
    !  - Calculate turbulent dry deposition velocity, vlc_trb.
    !  - Add vlc_trb to vlc_grv to give the total deposition velocity, vlc_dry.
    !------------------------------------------------------------------------------------
    call modal_aero_turb_drydep_velocity( moment, ncol, pcols, lchnk, aer_rmax,      &! in
                                          tair(:,pver), pmid(:,pver),                &! in
                                          radius_part(:,pver), density_part(:,pver), &! in
                                          sig_part(:,pver),                          &! in
                                          fricvel(:), ram1(:), vlc_grv(:,pver),      &! in
                                          vlc_trb(:), vlc_dry(:,pver)                )! out

  end subroutine modal_aero_depvel_part

  !=================================================================================
  ! Calculate some diagnostics for output. This does not affect time integration.
  !=================================================================================
  subroutine drydep_diags_for_1_tracer( lchnk, ncol, cnst_name_in, vlc_dry, vlc_trb, vlc_grv, sflx, dqdt_sed, pdel )

    use physconst,      only: gravit
    use cam_history,    only: outfld
    use cam_abortutils, only: endrun

    integer, intent(in) :: lchnk  ! chunk index
    integer, intent(in) :: ncol   ! # of active columns 

    character(len=*), intent(in) :: cnst_name_in  ! tracer name

    real(r8),intent(in) :: vlc_trb(pcols)       ! deposition velocity of turbulent dry deposition [m/s]
    real(r8),intent(in) :: vlc_grv(pcols,pver)  ! deposition velocity of gravitational settling [m/s]
    real(r8),intent(in) :: vlc_dry(pcols,pver)  ! deposition velocity of both mechanisms combined  [m/s]
    real(r8),intent(in) ::    sflx(pcols)       ! total deposition flux at the surface for one species [kg/m2/s] or [1/m2/s]

    real(r8),intent(in),optional :: dqdt_sed(pcols,pver)
    real(r8),intent(in),optional ::     pdel(pcols,pver)

    real(r8) :: dep_trb(pcols)       ! turbulent dry deposition portion of sflx [kg/m2/s] or [1/m2/s]
    real(r8) :: dep_grv(pcols)       ! gravitational settling   portion of slfx [kg/m2/s] or [1/m2/s]
    real(r8) :: tnd_trb(pcols)       ! diagnosed tendency corresponding to dep_trb [kg/kg/s] or [1/kg/s]
    real(r8) :: tnd_grv(pcols,pver)  ! diagnosed tendency corresponding to dep_trb [kg/kg/s] or [1/kg/s]
    integer :: ii

    !----------
    ! Fluxes
    !----------
    dep_trb(:) = 0._r8
    dep_grv(:) = 0._r8

    ! apportion dry deposition into turb and gravitational settling for tapes

    do ii=1,ncol
       if (vlc_dry(ii,pver) .ne. 0._r8) then
          dep_trb(ii)=sflx(ii)*vlc_trb(ii)/vlc_dry(ii,pver)
          dep_grv(ii)=sflx(ii)*vlc_grv(ii,pver)/vlc_dry(ii,pver)
       endif
    enddo

    ! send diagnostics to output

    call outfld( cnst_name_in//'DDF', sflx,     pcols, lchnk)
    call outfld( cnst_name_in//'TBF', dep_trb,  pcols, lchnk)
    call outfld( cnst_name_in//'GVF', dep_grv,  pcols, lchnk)

    !-------------
    ! Tendencies
    !-------------
    if ( present(dqdt_sed) ) then 

       if(.not.present(pdel)) call endrun("drydep_diags_for_1_tracer: please check input arguments")

       ! Diagnose tendency caused by turbulent dry deposition
       tnd_trb(1:ncol) = -dep_trb(1:ncol)*gravit/pdel(1:ncol,pver)

       ! Diagnose tendency caused by gravitational settling
       tnd_grv(1:ncol,:)    = dqdt_sed(1:ncol,:)
       tnd_grv(1:ncol,pver) = dqdt_sed(1:ncol,pver) - tnd_trb(1:ncol)

       ! Send to output
       call outfld( cnst_name_in//'DTQ_TB', tnd_trb,  pcols, lchnk)
       call outfld( cnst_name_in//'DTQ_GV', tnd_grv,  pcols, lchnk)
       call outfld( cnst_name_in//'DTQ',    dqdt_sed, pcols, lchnk)

     end if

  end subroutine drydep_diags_for_1_tracer

end module modal_aero_drydep_refac
