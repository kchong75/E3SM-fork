The branch `huiwanpnnl/gmd_2020_330_forc+cflx_202305` contains code modifications from 2023 that provide more options for the numerical coupling of the surface emissions, dry deposition, and turbulent mixing of
interstitial aerosols in EAMv1.

# Get the code

To get the code just for **viewing** the EAM part, use

```
git clone \
    -b huiwanpnnl/gmd_2020_330_forc+cflx_202305 \
    git@github.com:PAESCAL-SciDAC5/E3SM-fork.git \
    cflx_202305
```

If you want to **compile the branch and run simulations**, add `--recursive` to the clone command to also fetch all the submodules:

```
git clone --recursive \
    -b huiwanpnnl/gmd_2020_330_forc+cflx_202305 \
    git@github.com:PAESCAL-SciDAC5/E3SM-fork.git \
    cflx_202305
```

# View code changes on GitHub

Use the following URL in a web browser:

```
https://github.com/PAESCAL-SciDAC5/E3SM-fork/compare/v1_cflx_202305_starting_point...PAESCAL-SciDAC5:E3SM-fork:huiwanpnnl/gmd_2020_330_forc+cflx_202305?expand=1
```

# History of the branch

1. The branch was created off E3SM's main branch at `a2b742` in June 2018, several months after E3SMv1 was released.
1. Some coupling options for cloud processes were added until `a8de97` (December 9, 2021).
1. `clfx_cpl_opt = 1, 2, 3` were added and updated until `0dd149` (October 9, 2022). `clfx_cpl_opt = 1, 2` were integrated back to E3SM's main branch in March 2023 on top of v2 and in preparation for v3.
1. Starting in May 2023, more aerosol process coupling options were added, all controlled by the namelist variable `clfx_cpl_opt`. **The URL provided above shows only these code changes since May 2023.**

# Aerosol process coupling options

As of August 2023, two groups of options are available:

- `clfx_cpl_opt = 1, 2, 3, 4`. This group calculates gravitational settling and turbulent dry deposition of aerosols together, in `tphysac`, as in EAMv1 and v2. The 4 options differ in where/how the surface fluxes saved in `cam_in%cflx` are applied. The dry deposition code in use is the refactored one.

- `clfx_cpl_opt = 41, 42, 43, 44`. This group is the same as `clfx_cpl_opt = 4` in terms of where/how the surface fluxes saved in `cam_in%cflx` are applied, but unlike the first group, the turbulent dry deposition of interstitial aerosols is handled together with turbulent mixing, i.e., in `dropmixnuc`.

There is also `clfx_cpl_opt = 40` which was implemented for a sanity check. This option is numerically equivalent to `clfx_cpl_opt = 4` but uses the old dry deposition code, so that one can run two simulations using`clfx_cpl_opt = 4` and `40`, and then compare the global total energy/mass statistics in `atm.log.*` (i.e., the lines starting with "nstep, te") to verify BFB identify.

Detailed descriptions of the schemes can be found in [`note/cflx_cpl_opt_scheme_descriptions.tex`](note/cflx_cpl_opt_scheme_descriptions.tex)


# Overview of code changes

## BFB refactor of aerosol dry deposition code

The subroutine `aero_model_drydep` used in EAMv1 and v2, as well as some of the subroutines it calls, was refacotored. This was a necessary prepartion for separating graviational settling and turbulent dry deposition of interstitial aerosols.

- The old subroutine `aero_model_drydep` in `src/chemistry/modal_aero/aero_model.F90` is kept intact for reference.

- A new subroutine `aero_model_drydep_main`, which does the same numerical calculations as the old one and gives BFB identical results, can be found in `src/chemistry/modal_aero/modal_aero_drydep_refac.F90`.

- The following modules were added to `src/chemistry/modal_aero/` during the refactor:
  - `modal_aero_drydep_utils.F90`
  - `modal_aero_grav_setl.F90`
  - `modal_aero_turb_drydep.F90`
  - `mo_spitfire_transport.F90`

## Implementation of new schemes

In order to provide `clfx_cpl_opt = 4, 41, 42, 43, 44`, the following modules in `src/physics/cam/` were revised:

1. `physpkg.F90`: new options were added to `tphysac` and `tphysbc`.
1. `microp_aero.F90` and `ndrop.F90`: subroutine `dropmixnuc` can now also handle surface emissions and turbulent dry deposition of interstitial aerosols.
1. `sfc_cpl_opt.F90`: in subroutine `cflx_tend`, an option was added to save the aerosol mixing ratio tendencies to a separate array and zeroed out in `ptend`.

In addition, `addfld` and `add_default` calls for some new output variables have been added to subroutine `aero_model_init` in `src/chemistry/modal_aero/aero_model.F90`

## Code structure and diagnostics

See flow charts in [`notes/updated_drydep_output.drawio`](notes/updated_drydep_output.drawio). Note there are multiple pages/tabs in that file.
(**To open this type of files, please first download `draw.io` from [https://www.drawio.com](https://www.drawio.com)**.)

# Run scripts

See [scripts/](scripts)


