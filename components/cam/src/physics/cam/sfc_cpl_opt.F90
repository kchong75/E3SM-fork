module sfc_cpl_opt

  use shr_kind_mod,  only: r8=>shr_kind_r8
  use physconst,     only: gravit
  use ppgrid,        only: pcols, pver
  use constituents,  only: pcnst

  implicit none
  public

contains

  subroutine cflx_tend(state, cam_in, ptend, aero_cflx_tend)

    use physics_types,   only: physics_state, physics_ptend, physics_ptend_init
    use camsrfexch,      only: cam_in_t

    use phys_control,    only: phys_getopts
    use rad_constituents,only: rad_cnst_get_info, rad_cnst_get_mode_num_idx, rad_cnst_get_mam_mmr_idx

    implicit none

    type(physics_state), intent(in)     :: state                ! Physics state variables
    type(cam_in_t),      intent(in)     :: cam_in
    type(physics_ptend), intent(out)    :: ptend                ! Individual parameterization tendencies

    real(r8), intent(out), optional     :: aero_cflx_tend(pcols,pcnst)

    logical  :: lq(pcnst)
    integer  :: ncol, m
   !real(r8) :: tmp1(pcols)

    logical  :: prog_modal_aero
    integer  :: imode, ispec, icnst
    integer  :: nmodes,nspec

    ncol = state%ncol

   !----------
   !lq(:) = .TRUE.
   !call physics_ptend_init(ptend, state%psetcols, 'cflx_tend', lq=lq)

   !rztodt                 = 1._r8/ztodt
   !ptend%q(:ncol,:pver,:) = state%q(:ncol,:pver,:)
   !tmp1(:ncol)            = ztodt * gravit * state%rpdel(:ncol,pver)

   !do m = 2, pcnst
   !  ptend%q(:ncol,pver,m) = ptend%q(:ncol,pver,m) + tmp1(:ncol) * cam_in%cflx(:ncol,m)
   !enddo

   !ptend%q(:ncol,:pver,:) = (ptend%q(:ncol,:pver,:) - state%q(:ncol,:pver,:)) * rztodt
   !----------

    !----------
    lq(:) = .TRUE.
    call physics_ptend_init(ptend, state%psetcols, 'cflx_tend', lq=lq)
    
    do m = 2, pcnst
       ptend%q(:ncol,pver,m) = gravit * state%rpdel(:ncol,pver)* cam_in%cflx(:ncol,m)
    enddo
    !----------

    !------------------------------------------------------------------
    ! Move cflx-induced tendencies of aerosols to a different array
    ! so that they can be passed to dropmixnuc instead of being applied
    ! by the next call of physics_update.
    !------------------------------------------------------------------
    if (PRESENT(aero_cflx_tend)) then

       call phys_getopts( prog_modal_aero_out=prog_modal_aero )

       if (prog_modal_aero) then

           ! Get the number of aerosol modes
           call rad_cnst_get_info(0, nmodes=nmodes)

           ! Loop over all modes
           do imode = 1, nmodes

              ! Number mixing ratio of the mode
              call rad_cnst_get_mode_num_idx(imode, icnst)

              aero_cflx_tend(:ncol,icnst) = ptend%q(:ncol,pver,icnst)
              ptend%q(:ncol,:,icnst) = 0._r8
              lq(icnst)=.false.

              ! All mass mixing ratios of the mode

              call rad_cnst_get_info(0, imode, nspec=nspec)  ! # of species in the mode
              do ispec = 1, nspec
                 call rad_cnst_get_mam_mmr_idx(imode, ispec, icnst)
                 aero_cflx_tend(:ncol,icnst) = ptend%q(:ncol,pver,icnst)
                 ptend%q(:ncol,:,icnst) = 0._r8
                 lq(icnst)=.false.
              end do ! species in a mode

           end do ! modes

      end if ! if prog_modal_aero
    end if   ! if PRESENT(aero_cflx_tend)

  end subroutine cflx_tend

end module sfc_cpl_opt
