#!/bin/csh
date

set echo verbose

set cflx_cpl_opt  = 2 

set fetch_code    = 0   # 0 = No, >0 = Yes
set compile_model = 0   # 0 = No, >0 = Yes
set run_model     = 1   # 0 = No, >0 = Yes

####################################################################
# Fetch code
####################################################################
setenv CODEROOT $HOME/codes
setenv CCSMTAG  gmd_2020_330_cflx_202305
setenv CCSMROOT ${CODEROOT}/${CCSMTAG}

setenv BRANCH huiwanpnnl/gmd_2020_330_forc+cflx_202305
setenv REPO git@github.com:PAESCAL-SciDAC5/E3SM-fork.git

if ($fetch_code > 0) then

   cd ${CODEROOT}
   git clone -b $BRANCH --recursive $REPO $CCSMTAG

endif

####################################################################
# Machine, compset, PE layout etc.
####################################################################

setenv COMPSET F20TRC5-CMIP6 
setenv MRES ne30
setenv RESOLUTION ${MRES}_${MRES}
setenv MACH      compy
setenv PTMP      /compyfs/$user/scidac_scratch/2023_cflx/2023code

setenv ntasks 400
setenv nthrds 1

setenv CASE     cflx_cpl_opt${cflx_cpl_opt}_bgt_${COMPSET}

setenv CASEROOT  $PTMP/$CASE/case
setenv RUNDIR    $PTMP/$CASE/run

#
# RUNDIR: $MEMBERWORK/$PROJECT/$CASE/run
# EXEROOT: $CESMSCRATCHROOT/$CASE/bld
# CESMSCRATCHROOT: $HOME/acme_scratch/$PROJECT
#
####################################################################
# Compile model
####################################################################
if ($compile_model > 0 || $run_model > 0) then

   rm -rf $CASEROOT
   cd  $CCSMROOT/cime/scripts

###usage: create_newcase [-h] [-d] [-v] [-s] --case CASE --compset COMPSET --res
###                      RES [--mach MACH] [--compiler COMPILER] [--ninst NINST]
###                      [--mpilib MPILIB] [--project PROJECT]
###                      [--pecount PECOUNT] [--mach-dir MACH_DIR]
###                      [--user-mods-dir USER_MODS_DIR] [--user-compset]
###                      [--pesfile PESFILE] [--user-grid] [--gridfile GRIDFILE]
###                      [--srcroot SRCROOT] [--test] [--walltime WALLTIME]
###                      [-q QUEUE]


   ./create_newcase --case $CASEROOT --project esmd --mach $MACH \
                    --res $RESOLUTION --compset $COMPSET

#====================================================================
# set up case
#====================================================================

   ###./create_newcase -list grids

   cd $CASEROOT

   ./xmlchange -file env_run.xml      -id RUNDIR     -val $RUNDIR
   ./xmlchange -file env_build.xml    -id EXEROOT    -val $PTMP/cflx_BLD1_$CCSMTAG

   ./xmlchange -file env_mach_pes.xml -id NTASKS_ATM -val $ntasks
   ./xmlchange -file env_mach_pes.xml -id NTHRDS_ATM -val $nthrds
   ./xmlchange -file env_mach_pes.xml -id ROOTPE_ATM -val '0'

   ./xmlchange -file env_mach_pes.xml -id NTASKS_LND -val $ntasks
   ./xmlchange -file env_mach_pes.xml -id NTHRDS_LND -val $nthrds
   ./xmlchange -file env_mach_pes.xml -id ROOTPE_LND -val '0'

   ./xmlchange -file env_mach_pes.xml -id NTASKS_ROF -val $ntasks
   ./xmlchange -file env_mach_pes.xml -id NTHRDS_ROF -val $nthrds
   ./xmlchange -file env_mach_pes.xml -id ROOTPE_ROF -val '0'

   ./xmlchange -file env_mach_pes.xml -id NTASKS_ICE -val $ntasks
   ./xmlchange -file env_mach_pes.xml -id NTHRDS_ICE -val $nthrds
   ./xmlchange -file env_mach_pes.xml -id ROOTPE_ICE -val '0'

   ./xmlchange -file env_mach_pes.xml -id NTASKS_OCN -val $ntasks
   ./xmlchange -file env_mach_pes.xml -id NTHRDS_OCN -val $nthrds
   ./xmlchange -file env_mach_pes.xml -id ROOTPE_OCN -val '0'

   ./xmlchange -file env_mach_pes.xml -id NTASKS_GLC -val $ntasks
   ./xmlchange -file env_mach_pes.xml -id NTHRDS_GLC -val $nthrds
   ./xmlchange -file env_mach_pes.xml -id ROOTPE_GLC -val '0'

   ./xmlchange -file env_mach_pes.xml -id NTASKS_WAV -val $ntasks
   ./xmlchange -file env_mach_pes.xml -id NTHRDS_WAV -val $nthrds
   ./xmlchange -file env_mach_pes.xml -id ROOTPE_WAV -val '0'

   ./xmlchange -file env_mach_pes.xml -id NTASKS_CPL -val $ntasks
   ./xmlchange -file env_mach_pes.xml -id NTHRDS_CPL -val $nthrds
   ./xmlchange -file env_mach_pes.xml -id ROOTPE_CPL -val '0'

   ./case.setup

#====================================================================
# Compile 
#====================================================================
   cd $CASEROOT

#  ./xmlchange -file env_build.xml -id DEBUG       -val 'TRUE'
   ./xmlchange -file env_build.xml -id PIO_VERSION -val '1'

   if ($compile_model > 0) then # Build the model

      ./case.build

   else # mark as already built

      ./xmlchange -file env_build.xml -id BUILD_COMPLETE  -val 'TRUE'

   endif

endif

#####################################################################
# Conduct simulation
#####################################################################
if ($run_model > 0) then
#------------------
## set environment
#------------------

cd $CASEROOT

./xmlchange  -file env_run.xml  -id  RUN_STARTDATE      -val '2009-10-01'
./xmlchange  -file env_run.xml  -id  RESUBMIT           -val '0'
./xmlchange  -file env_run.xml  -id  CONTINUE_RUN       -val 'FALSE'
./xmlchange  -file env_run.xml  -id  STOP_N             -val '3'
./xmlchange  -file env_run.xml  -id  STOP_OPTION        -val 'ndays'
./xmlchange  -file env_run.xml  -id  REST_N             -val '3'
./xmlchange  -file env_run.xml  -id  REST_OPTION        -val 'ndays'
./xmlchange  -file env_run.xml  -id  DOUT_S             -val 'FALSE'

./xmlchange  -file env_run.xml  -id  PIO_TYPENAME       -val 'netcdf'

./xmlchange  -file env_batch.xml -id JOB_WALLCLOCK_TIME -val '1:00:00'
#./xmlchange JOB_QUEUE=short --force

#./xmlchange  -file env_batch.xml -id JOB_QUEUE          -val 'short'

cat <<EOF >! user_nl_cam
&camexp
!
cflx_cpl_opt = ${cflx_cpl_opt}
!
!...................
! conditional diag
!...................
metric_name = 'ALL',
!
qoi_chkpt = 'CHEM','CFLX1','AERDRYRM','PBCINI','CFLX2'
            'STCLD','AERWETRM',
!
qoi_name = 'dst_a1','dst_a1','dst_a3','dst_a3'
qoi_nver =  72,      72,      72,      72
qoi_x_dp =   0,       2,       0,       2
!
l_output_state = .true.
l_output_incrm = .true.
!
!
!.......................................................
! history files
!.......................................................
 hist_tape_with_all_output = 1,

 nhtfrq          =  -24,
 mfilt           =  1,
 avgflag_pertape = 'A'

!----
!hist_tape_with_all_output = 1,2,3
!
!fincl2lonlat = '0e:360e_5n:55n'
!fincl3lonlat = '0e:360e_5n:55n'
!
!nhtfrq          =  0,  0 ,  -6,
!mfilt           =  1,  1   120,
!avgflag_pertape = 'A', 'A','I',
!----

 history_amwg        = .false.
 history_aero_optics = .false.
 history_aerosol     = .true.
 history_verbose     = .true.
!
!...................
! change init data
!...................
ncdata = '/compyfs/sunj695/csmruns/init/constance_FC5AV1C-04P2_ne30_ne30_E3SMv1.0_CLIMO_PD_TUNE.cam.i.0000-10-01-00000.nc'
/
EOF

cat <<EOB >! user_nl_clm
&clm_inparm
 finidat = '/compyfs/sunj695/csmruns/init/constance_FC5AV1C-04P2_ne30_ne30_E3SMv1.0_CLIMO_PD_TUNE.clm2.r.0000-10-01-00000.nc'
 check_finidat_fsurdat_consistency = .false.
/
EOB

./case.submit

endif
