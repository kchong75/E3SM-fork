!===============================================================================
! Calculations of 
!  - changes in interstitial aerosol mixing ratios caused by
!    gravitational settling and turbulent dry deposition of aerosol particles,
!  - changes in cloud-borne aerosol mixing ratios caused by 
!    gravitational settling and turbulent dry deposition of cloud droplets.
!===============================================================================
module modal_aero_turb_drydep

  use shr_kind_mod,   only: r8 => shr_kind_r8
  use ppgrid,         only: pcols
  use physconst,      only: gravit

  implicit none
  public

contains
  
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
