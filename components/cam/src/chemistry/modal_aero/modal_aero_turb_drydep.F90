!===============================================================================
! Calculations of 
!  - changes in interstitial aerosol mixing ratios caused by
!    gravitational settling and turbulent dry deposition of aerosol particles,
!  - changes in cloud-borne aerosol mixing ratios caused by 
!    gravitational settling and turbulent dry deposition of cloud droplets.
!===============================================================================
module modal_aero_turb_drydep

  use shr_kind_mod,   only: r8 => shr_kind_r8
  use ppgrid,         only: pcols, pver
  use physconst,      only: gravit, rair
  use perf_mod,       only: t_startf, t_stopf
  use cam_history,    only: outfld

  implicit none
  private

  public :: get_gridcell_ram1_fricvel
  public :: interstitial_aero_turb_dep_velocity
  public :: modal_aero_turb_drydep_velocity

contains

  !---------------------------------------------------------------------------------
  subroutine get_gridcell_ram1_fricvel( state, cam_in, ram1, fricvel )

    use physics_types,          only: physics_state
    use camsrfexch,             only: cam_in_t
    use clubb_intr,             only: calc_ustar_obklen => clubb_surface

    type(physics_state),    intent(in) :: state
    type(cam_in_t),         intent(in) :: cam_in

    real(r8),intent(out) ::    ram1(pcols) ! aerodynamical resistance used in the calculaiton of turbulent dry deposition velocity [s/m]
    real(r8),intent(out) :: fricvel(pcols) ! friction velocity used in the calculaiton of turbulent dry deposition velocity [m/s]

    real(r8) ::  ustar(pcols)      ! surface friction velocity
    real(r8) ::  obklen(pcols)     ! Obukhov length

    call calc_ustar_obklen( state, cam_in, ustar, obklen )
    call calc_ram( state%ncol,                                              &! in
                   cam_in%landfrac, cam_in%icefrac, cam_in%ocnfrac,         &! in
                   obklen,          ustar,                                  &! in; calculated above in tphysac
                   state%t(:,pver), state%pmid(:,pver), state%pdel(:,pver), &! in; note: bottom level only
                   cam_in%ram1,     cam_in%fv,                              &! in
                   ram1,            fricvel                               )  ! out

  end subroutine get_gridcell_ram1_fricvel

  !----------------------------------------------------------------------------------------------
  ! !DESCRIPTION: 
  !  
  ! Determine a single aerodynamic resistance value for the grid cell.
  !  - If the grid cell has non-zero land fraction, use the value from the coupler (ram1_in);
  !  - Otherwise (for oceans and sea ice), calculate a value from Seinfeld and Pandis
  !    (2006, 2nd edition, ISBN-13: 978-0471720188. p.907, Eq. 19.13-19.14)
  !  
  ! Determine a single friction velocity for the grid cell.
  !  - If the grid cell has non-zero land fraction, use the value from the coupler (fv_in);
  !  - Otherwise (for oceans and sea ice), use the value from the atmosphere model (ustar).
  !
  ! Author: Natalie Mahowald
  ! Code refactor: Hui Wan, 2023
  !----------------------------------------------------------------------------------------------
  subroutine calc_ram(ncol,landfrac,icefrac,ocnfrac, &
                     obklen,ustar,                  &
                     tair, pmid, pdel,              &
                     ram1_in, fv_in,                &
                     ram1_out,fv_out )

     implicit none
     integer, intent(in) :: ncol

     real(r8), intent(in) :: landfrac(pcols)  ! land fraction [unitless]
     real(r8), intent(in) :: icefrac(pcols)   ! ice fraction  [unitless]
     real(r8), intent(in) :: ocnfrac(pcols)   ! ocean fraction [unitless]

     real(r8), intent(in) :: obklen(pcols)    ! Obukhov length [m]
     real(r8), intent(in) :: ustar(pcols)     ! Surface friction velocity [m/s]

     real(r8), intent(in) :: tair(pcols)      ! air temperature [K]
     real(r8), intent(in) :: pmid(pcols)      ! air pressure [Pa]
     real(r8), intent(in) :: pdel(pcols)      ! layer pressure thickness [Pa]

     real(r8),intent(in) :: ram1_in(pcols)     ! aerodynamical resistance [s/m]
     real(r8),intent(in) :: fv_in(pcols)       ! sfc friction velocity from land model [m/s]

     real(r8),intent(out) :: ram1_out(pcols)   ! aerodynamical resistance [s/m]
     real(r8),intent(out) :: fv_out(pcols)     ! sfc friction velocity of the current grid cell [m/s]

     real(r8), parameter :: zzocen = 0.0001_r8   ! Ocean aerodynamic roughness length
     real(r8), parameter :: zzsice = 0.0400_r8   ! Sea ice aerodynamic roughness length
     real(r8), parameter :: xkar   = 0.4_r8      ! Von Karman constant

     ! local variables

     real(r8) :: zz,psi,psi0,nu,nu0,temp
     integer :: ii

     real(r8), parameter :: lndfrc_threshold = 0.000000001_r8  ! fraction, unitless

     !---------------------------------------------------------------------------
     ! Friction velocity:
     !  - If the grid cell has a land fraction larger than a threshold (~zero),
     !    then use fv_in (which comes from the coupler).
     !  - Otherwise, use the ustar calculated in the atmosphere.
     !---------------------------------------------------------------------------
     where( landfrac(:ncol) > lndfrc_threshold )
       fv_out(:ncol) = fv_in(:ncol)
     elsewhere
       fv_out(:ncol) = ustar(:ncol)
     endwhere

     ! fvitt -- fv == 0 causes a floating point exception in
     ! dry dep of sea salts and dust

     where ( fv_out(:ncol) == 0._r8 )
        fv_out(:ncol) = 1.e-12_r8
     endwhere

     !-------------------------------------------------------------------
     ! Aerodynamic resistence
     !-------------------------------------------------------------------
     do ii=1,ncol

        ! If the grid cell has a land fraction larger than a threshold (~zero),
        ! simply use ram1_in (which comes from the coupler)

        if (landfrac(ii) > lndfrc_threshold) then
           ram1_out(ii)=ram1_in(ii)

        else
        ! If the grid cell has a land fraction smaller than the threshold,
        ! calculate aerodynamic resistence


           ! calculate psi, psi0, temp

           zz=pdel(ii)*rair*tair(ii)/pmid(ii)/gravit/2.0_r8   !use half the layer height like Ganzefeld and Lelieveld, 1995
           if(obklen(ii).eq.0) then
              psi=0._r8
              psi0=0._r8
           else
              psi=min(max(zz/obklen(ii),-1.0_r8),1.0_r8)
              psi0=min(max(zzocen/obklen(ii),-1.0_r8),1.0_r8)
           endif

           temp=zz/zzocen

           ! special treatment for ice-dominant cells

           if(icefrac(ii) > 0.5_r8) then 
              if(obklen(ii).gt.0) then 
                 psi0=min(max(zzsice/obklen(ii),-1.0_r8),1.0_r8)
              else
                 psi0=0.0_r8
              endif
              temp=zz/zzsice
           endif

           ! calculate aerodynamic resistence

           if(psi> 0._r8) then
              ram1_out(ii) =1/xkar/ustar(ii)*(log(temp)+4.7_r8*(psi-psi0))
           else
              nu=(1.00_r8-15.000_r8*psi)**(.25_r8)
              nu0=(1.000_r8-15.000_r8*psi0)**(.25_r8)
              if(ustar(ii).ne.0._r8) then
                 ram1_out(ii) =1/xkar/ustar(ii)*(log(temp) &
                      +log(((nu0**2+1.00_r8)*(nu0+1.0_r8)**2)/((nu**2+1.0_r8)*(nu+1.00_r8)**2)) &
                      +2.0_r8*(atan(nu)-atan(nu0)))
              else
                 ram1_out(ii) =0._r8
              endif
           endif

        endif  ! if grid cell has land fraction > threshold or not

     enddo ! loop over grid columns

  end subroutine calc_ram

  !========================================================================================
  !========================================================================================
  subroutine interstitial_aero_turb_dep_velocity( state, pbuf, ram1, fricvel, vlc_grv, vlc_trb )

    use physics_types,    only: physics_state
    use physics_buffer,   only: physics_buffer_desc, pbuf_get_field
    use modal_aero_data,  only: ntot_amode
    use modal_aero_data,  only: alnsg_amode, sigmag_amode
    use aero_model,       only: dgnumwet_idx, wetdens_ap_idx, nmodes

    use modal_aero_drydep_utils, only: aer_rmax

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
       call modal_aero_turb_drydep_velocity( imnt, ncol, pcols, lchnk, aer_rmax,    &! in
                                             tair, pmid,                            &! in
                                             rad_aer, dens_aer, sg_aer,             &! in
                                             fricvel, ram1,                         &! in
                                             vlc_grv(:,jvlc,imode), &! in
                                             vlc_trb(:,jvlc,imode), &! out
                                             vlc_dry(:)             )! out

       jvlc = 2  ; imnt = 3  ! interstitial aerosol volume/mass
       call modal_aero_turb_drydep_velocity( imnt, ncol, pcols, lchnk, aer_rmax,    &! in
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
  
  !------------------------------------------------------------------------------------
  ! Calculate particle velocity of turbulent dry deposition.
  ! This process is assumed to only occur in the lowest model layer.
  !------------------------------------------------------------------------------------
  subroutine modal_aero_turb_drydep_velocity( moment, ncol, pcols, lchnk, radius_max,  &! in
                                              tair, pmid,                              &! in
                                              radius_part, density_part, sig_part,     &! in
                                              fricvel, ram1, vlc_grv,                  &! in
                                              vlc_trb, vlc_dry                         )! out

    use mo_drydep,     only: n_land_type, fraction_landuse
    use modal_aero_drydep_utils, only: air_dynamic_viscosity, air_kinematic_viscosity
    use modal_aero_drydep_utils, only: radius_for_moment, schmidt_number

    implicit none

    integer,  intent(in) :: moment           ! moment of size distribution (0 for number, 2 for surface area, 3 for volume)
    integer,  intent(in) :: ncol             ! # of grid columns to do calculations for
    integer,  intent(in) :: pcols            ! dimension size (# of grid columns)
    integer,  intent(in) :: lchnk            ! index of grid chunk

    real(r8), intent(in) :: radius_max       ! upper bound of radius used for the calculation of deposition velocity [m]

    real(r8), intent(in) :: tair(pcols)      ! air temperature [K]
    real(r8), intent(in) :: pmid(pcols)      ! air pressure [Pa]

    real(r8), intent(in) ::  radius_part(pcols)   ! mean (volume or number) particle radius [m]
    real(r8), intent(in) :: density_part(pcols)   ! density of particle material [kg/m3]
    real(r8), intent(in) ::     sig_part(pcols)   ! geometric standard deviation of particle size distribution

    real(r8), intent(in) :: fricvel(pcols)   ! friction velocity [m/s]
    real(r8), intent(in) :: ram1(pcols)      ! aerodynamical resistance [s/m]

    real(r8), intent(in)  :: vlc_grv(pcols)  ! gravitational deposition velocity [m/s]
    real(r8), intent(out) :: vlc_trb(pcols)  ! turbulent dry deposition velocity [m/s]
    real(r8), intent(out) :: vlc_dry(pcols)  ! total     dry deposition velocity [m/s]

    !------------------------------------------------------------------------
    ! Local Variables

    integer  :: ii            ! grid column index
    real(r8) :: vsc_dyn_atm   ! Dynamic viscosity of air [kg m-1 s-1]
    real(r8) :: vsc_knm_atm   ! Kinematic viscosity of atmosphere [m2 s-1]
    real(r8) :: slp_crc       ! Slip correction factor [unitless]
    real(r8) :: radius_moment ! median radius for moment [m]

    real(r8) :: shm_nbr       ! Schmidt number [unitless]
    real(r8) :: stk_nbr       ! Stokes number [unitless]
    real(r8) :: rss_trb       ! Resistance to turbulent deposition [s m-1]
    real(r8) :: rss_lmn       ! Quasi-laminar layer resistance [s m-1]
    real(r8) :: brownian      ! collection efficiency for Browning diffusion
    real(r8) :: impaction     ! collection efficiency for impaction
    real(r8) :: interception  ! collection efficiency for interception
    real(r8) :: stickfrac     ! fraction of particles sticking to surface

    integer  :: lt              ! land type index
    real(r8) :: lnd_frc         ! land type fraction [unitless]
    real(r8) :: vlc_trb_ontype  ! turbulent dry dep. velocity on single land type [m/s]
    real(r8) :: vlc_trb_wgtsum  ! turbulent dry dep. velocity averaged over land types [m/s]
    real(r8) :: vlc_dry_wgtsum  ! total     dry dep. velocity averaged over land types [m/s]

    ! constants

    real(r8),parameter :: eps0 = 3.0_r8 ! empirical parameter $\varepsilon_0$ in Eq. (5) of Zhang L. et al. (2001)
    real(r8),parameter :: beta = 2.0_r8 ! empirical parameter $\beta$ in Eq. (7c) of Zhang L. et al. (2001)
    real(r8),parameter :: stickfrac_lowerbnd  = 1.0e-10_r8  ! lower bound of stick fraction

    real(r8) gamma(11)      ! exponent of schmidt number
    data gamma/0.56e+00_r8,  0.54e+00_r8,  0.54e+00_r8,  0.56e+00_r8,  0.56e+00_r8, &        
               0.56e+00_r8,  0.50e+00_r8,  0.54e+00_r8,  0.54e+00_r8,  0.54e+00_r8, &
               0.54e+00_r8/
    save gamma

    real(r8) alpha(11)      ! parameter for impaction
    data alpha/1.50e+00_r8,   1.20e+00_r8,  1.20e+00_r8,  0.80e+00_r8,  1.00e+00_r8, &
               0.80e+00_r8, 100.00e+00_r8, 50.00e+00_r8,  2.00e+00_r8,  1.20e+00_r8, &
              50.00e+00_r8/
    save alpha

    real(r8) radius_collector(11) ! radius (m) of surface collectors
    data radius_collector/10.00e-03_r8,  3.50e-03_r8,  3.50e-03_r8,  5.10e-03_r8,  2.00e-03_r8, &
                           5.00e-03_r8, -1.00e+00_r8, -1.00e+00_r8, 10.00e-03_r8,  3.50e-03_r8, &
                          -1.00e+00_r8/
    save radius_collector

    integer  :: iwet(11) ! flag for wet surface = 1, otherwise = -1
    data iwet/-1,  -1,   -1,   -1,   -1,  &
              -1,   1,   -1,    1,   -1,  &
              -1/
    save iwet

    !---------
    do ii=1,ncol

       ! Calculate size-INdependent thermokinetic properties of the air

       vsc_dyn_atm = air_dynamic_viscosity( tair(ii) )
       vsc_knm_atm = air_kinematic_viscosity( tair(ii), pmid(ii) )

       ! Calculate the mean radius and Schmidt number of the moment

       radius_moment = radius_for_moment( moment,sig_part(ii),radius_part(ii),radius_max )
       shm_nbr = schmidt_number( tair(ii), pmid(ii), radius_moment, vsc_dyn_atm, vsc_knm_atm )

       ! Initialize deposition velocities averages over different land surface types

       vlc_trb_wgtsum = 0._r8
       vlc_dry_wgtsum = 0._r8

       ! Loop over different land surface types. Calculate deposition velocities of 
       ! those different surface types. The overall deposition velocity of a grid cell 
       ! is the area-weighted average of those land-type-specific velocities.

       do lt = 1,n_land_type

          lnd_frc = fraction_landuse(ii,lt,lchnk)

          if ( lnd_frc /= 0._r8 ) then

             !----------------------------------------------------------------------
             ! Collection efficiency of deposition mechanism 1 - Brownian diffusion
             !----------------------------------------------------------------------
             brownian = shm_nbr**(-gamma(lt))

             !----------------------------------------------------------------------
             ! Collection efficiency of deposition mechanism 2 - interception 
             !----------------------------------------------------------------------
             if (radius_collector(lt) > 0.0_r8) then ! vegetated surface
                interception = 2.0_r8*(radius_moment/radius_collector(lt))**2.0_r8
             else ! non-vegetated surface
                interception = 0.0_r8
             endif

             !----------------------------------------------------------------------
             ! Collection efficiency of deposition mechanism 3 - impaction 
             !----------------------------------------------------------------------
             if (radius_collector(lt) > 0.0_r8) then ! vegetated surface
                stk_nbr = vlc_grv(ii) * fricvel(ii) / (gravit*radius_collector(lt))
             else ! non-vegetated surface
                stk_nbr = vlc_grv(ii) * fricvel(ii) * fricvel(ii) / (gravit*vsc_knm_atm)  ! SeP97 p.965
             endif

             impaction = (stk_nbr/(alpha(lt)+stk_nbr))**beta   ! Eq. (7c) of Zhang L. et al.  (2001)

             !-----------------------------------------------------
             ! Stick fraction, Eq. (10) of Zhang L. et al.  (2001)
             !-----------------------------------------------------
             if (iwet(lt) > 0) then
                stickfrac = 1.0_r8
             else
                stickfrac = max( stickfrac_lowerbnd, exp(-sqrt(stk_nbr)) )
             endif

             !----------------------------------------------------------------------------------
             ! Using the numbers calculated above, compute the quasi-laminar layer resistance
             ! following Zhang L. et al. (2001), Eq. (5)
             !----------------------------------------------------------------------------------
             rss_lmn = 1.0_r8 / (eps0 * fricvel(ii) * stickfrac * (brownian+interception+impaction))

             !--------------------------------------------------------------------
             ! Total resistence and deposition velocity of turbulent deposition,
             ! see Eq. (21) of Zender et al. (2003)
             !--------------------------------------------------------------------
             rss_trb = ram1(ii) + rss_lmn + ram1(ii)*rss_lmn*vlc_grv(ii)
             vlc_trb_ontype = 1.0_r8 / rss_trb

             !--------------------------------------------------------------------------------
             ! Contributions to the single-value bulk deposition velocities of the grid cell
             !--------------------------------------------------------------------------------
             vlc_trb_wgtsum = vlc_trb_wgtsum + lnd_frc*( vlc_trb_ontype )
             vlc_dry_wgtsum = vlc_dry_wgtsum + lnd_frc*( vlc_trb_ontype + vlc_grv(ii) )
          endif

       enddo  ! lt=1,n_land_type

       vlc_trb(ii) = vlc_trb_wgtsum
       vlc_dry(ii) = vlc_dry_wgtsum

    enddo ! ii=1,ncol
    !-------

  end subroutine modal_aero_turb_drydep_velocity

end module modal_aero_turb_drydep
