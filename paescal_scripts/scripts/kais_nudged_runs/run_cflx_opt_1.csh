#!/bin/csh
date

set echo verbose

set compile_model = 1   # 0 = No, >0 = Yes
set run_model     = 1   # 0 = No, >0 = Yes
set walltime      = '7:59:00'

####################################################################
# Fetch code
####################################################################
setenv E3SMTAG  gmd_2020_330_cflx
setenv E3SMROOT $HOME/compy/model/${E3SMTAG}

####################################################################
# Machine, compset, PE layout etc.
####################################################################

setenv COMPSET FC5AV1C-L 
setenv MRES ne30
setenv RESOLUTION ${MRES}_${MRES}
setenv MACH      compy
setenv PTMP      /compyfs/$user/e3sm_scratch/

setenv NTASKS 1200
setenv NTHRDS 1

setenv CASE     ${MACH}_${COMPSET}_${MRES}_ndg_cflx_cpl_opt1

setenv CASEROOT  ${E3SMROOT}/cases/$CASE
setenv RUNDIR    $PTMP/$CASE/run

#
# set RUNDIR     = /compyfs/zhan524/e3sm_scratch/$CASE
#

####################################################################
# Compile model
####################################################################
if ($compile_model > 0 || $run_model > 0) then

   rm -rf $CASEROOT
   cd  $E3SMROOT/cime/scripts

###usage: create_newcase [-h] [-d] [-v] [-s] --case CASE --compset COMPSET --res
###                      RES [--mach MACH] [--compiler COMPILER] [--ninst NINST]
###                      [--mpilib MPILIB] [--project PROJECT]
###                      [--pecount PECOUNT] [--mach-dir MACH_DIR]
###                      [--user-mods-dir USER_MODS_DIR] [--user-compset]
###                      [--pesfile PESFILE] [--user-grid] [--gridfile GRIDFILE]
###                      [--srcroot SRCROOT] [--test] [--walltime WALLTIME]
###                      [-q QUEUE]


   ./create_newcase --case $CASEROOT --project e3sm --mach $MACH \
                    --res $RESOLUTION --compset $COMPSET

####################################################################
# set up case
####################################################################

#./create_newcase -list grids

   cd $CASEROOT

#  ./xmlchange RUNDIR=$RUNDIR
   ./xmlchange EXEROOT=$PTMP/$CASE/bld/

   ./xmlchange NTASKS_ATM=$NTASKS
   ./xmlchange NTHRDS_ATM=$NTHRDS
   ./xmlchange ROOTPE_ATM='0'

   ./xmlchange NTASKS_LND=$NTASKS
   ./xmlchange NTHRDS_LND=$NTHRDS
   ./xmlchange ROOTPE_LND='0'

   ./xmlchange NTASKS_ROF=$NTASKS
   ./xmlchange NTHRDS_ROF=$NTHRDS
   ./xmlchange ROOTPE_ROF='0'

   ./xmlchange NTASKS_ICE=$NTASKS
   ./xmlchange NTHRDS_ICE=$NTHRDS
   ./xmlchange ROOTPE_ICE='0'

   ./xmlchange NTASKS_OCN=$NTASKS
   ./xmlchange NTHRDS_OCN=$NTHRDS
   ./xmlchange ROOTPE_OCN='0'

   ./xmlchange NTASKS_GLC=$NTASKS
   ./xmlchange NTHRDS_GLC=$NTHRDS
   ./xmlchange ROOTPE_GLC='0'

   ./xmlchange NTASKS_WAV=$NTASKS
   ./xmlchange NTHRDS_WAV=$NTHRDS
   ./xmlchange ROOTPE_WAV='0'

   ./xmlchange NTASKS_CPL=$NTASKS
   ./xmlchange NTHRDS_CPL=$NTHRDS
   ./xmlchange ROOTPE_CPL='0'

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

./xmlchange  JOB_QUEUE=slurm --force
./xmlchange  JOB_WALLCLOCK_TIME=$walltime
./xmlchange  RUN_STARTDATE='2009-10-01'
./xmlchange  RESUBMIT=4
./xmlchange  STOP_N='3'
./xmlchange  STOP_OPTION='nmonths'
./xmlchange  REST_N='3'
./xmlchange  REST_OPTION='nmonths'
#./xmlchange CONTINUE_RUN='TRUE'
#./xmlchange DOUT_S='FALSE'
#./xmlchange DOUT_L_MS='FALSE'

##./xmlchange SSTICE_DATA_FILENAME='/lustre/atlas1/cli112/scratch/kaizhang/input/sst_HadOIBl_bc_1x1_clim_pi_plus4K.nc'


#####################################################################
# Namelist changes  
#####################################################################
cat <<EOF >! user_nl_cam
&camexp
inithist = 'MONTHLY'
inithist_all = .true.
!.......................................................
! history files
!.......................................................
 nhtfrq          =  0, !!  -24
 mfilt           =  1, !!   31
 avgflag_pertape = 'A' !!, 'A'
!
 rad_diag_1 = 'A:Q:H2O', 'N:O2:O2', 'N:CO2:CO2', 
 	     'A:O3:O3', 'N:N2O:N2O', 'N:CH4:CH4', 
 	     'N:CFC11:CFC11', 'N:CFC12:CFC12',
 docosp    = .true.,
 cosp_lite = .true.,
 cosp_ncolumns        = 10
 cosp_nradsteps       = 3
 cosp_lmisr_sim       = .true.
 cosp_lisccp_sim      = .true.
 cosp_lmodis_sim      = .true.
 cosp_llidar_sim      = .true.
 history_amwg         = .true.
 history_aero_optics  = .true.
 history_aerosol      = .true.
 history_clubb        = .true.
 history_budget       = .true.
 history_verbose      = .true.
 hist_hetfrz_classnuc = .true. 
 do_aerocom_ind3 = .true.
!...................
! change init data
!...................
ncdata = '/compyfs/zhan524/myinput/init/compy_FC5AV1C-04P2_ne30_E3SM_20190418_ctrl_1600p_EXP01.cam.i.0000-10-01-00000.nc'
!.......................................................
! nudging
!
! Default setup in E3SMv1 
!    Nudge_Tau         = -999
!    Nudge_Loc_PhysOut = .False.
!    Nudge_CurrentStep = .False.
!    Nudge_File_Ntime  = 1
!    Nudge_Method      = ‘Step’
!
! Setup for MERRA2 
!    Nudge_Path           = /qfs/projects/eagles/pma/merra2/ne30np4/'
!    Nudge_File_Template  = 'MERRA2_ne30np4_%y-%m-%d-%s.nc'
!    Nudge_Times_Per_Day  = 8  !! nudging input data frequency
!    Nudge_File_Ntime     = 1 
!.......................................................
Nudge_Model          = .True.
Nudge_Path           = '/compyfs/zhan524/TMP/ndata/eraint_ne30L72_se/'
Nudge_File_Template  = 'interim_se_%y%m%d00.nc'
Nudge_Times_Per_Day  = 4  !! nudging input data frequency
Model_Times_Per_Day  = 48 !! should not be larger than 48 if dtime = 1800s
Nudge_Uprof          = 2
Nudge_Ucoef          = 1.
Nudge_Vprof          = 2
Nudge_Vcoef          = 1.
Nudge_Tprof          = 0
Nudge_Tcoef          = 0.
Nudge_Qprof          = 0
Nudge_Qcoef          = 0.
Nudge_PSprof         = 0
Nudge_PScoef         = 0.
Nudge_Beg_Year       = 0001
Nudge_Beg_Month      = 1
Nudge_Beg_Day        = 1
Nudge_End_Year       = 9999
Nudge_End_Month      = 1
Nudge_End_Day        = 1
Nudge_Vwin_Lindex    = 0.
Nudge_Vwin_Hindex    = 70.
Nudge_Vwin_Ldelta    = 0.1
Nudge_Vwin_Hdelta    = 0.1
Nudge_Vwin_lo        = 0.
Nudge_Vwin_hi        = 1.
Nudge_Method         = 'Linear'
Nudge_Loc_PhysOut    = .True.
Nudge_Tau            = 6.        !! relaxation time scale, unit: 6h
Nudge_CurrentStep    = .False.
Nudge_File_Ntime     = 4
!.......................................................
! emissions 
!.......................................................
ext_frc_specifier              = 'SO2         -> /compyfs/inputdata/atm/cam/chem/trop_mozart_aero/emis/DECK_ne30/cmip6_mam4_so2_elev_1850-2014_c180205.nc',
        'SOAG        -> /compyfs/inputdata/atm/cam/chem/trop_mozart_aero/emis/DECK_ne30/cmip6_mam4_soag_elev_1850-2014_c180205.nc',
        'bc_a4       -> /compyfs/inputdata/atm/cam/chem/trop_mozart_aero/emis/DECK_ne30/cmip6_mam4_bc_a4_elev_1850-2014_c180205.nc',
        'num_a1      -> /compyfs/inputdata/atm/cam/chem/trop_mozart_aero/emis/DECK_ne30/cmip6_mam4_num_a1_elev_1850-2014_c180205.nc',
        'num_a2      -> /compyfs/inputdata/atm/cam/chem/trop_mozart_aero/emis/DECK_ne30/cmip6_mam4_num_a2_elev_1850-2014_c180205.nc',
        'num_a4      -> /compyfs/inputdata/atm/cam/chem/trop_mozart_aero/emis/DECK_ne30/cmip6_mam4_num_a4_elev_1850-2014_c180205.nc',
        'pom_a4      -> /compyfs/inputdata/atm/cam/chem/trop_mozart_aero/emis/DECK_ne30/cmip6_mam4_pom_a4_elev_1850-2014_c180205.nc',
        'so4_a1      -> /compyfs/inputdata/atm/cam/chem/trop_mozart_aero/emis/DECK_ne30/cmip6_mam4_so4_a1_elev_1850-2014_c180205.nc',
        'so4_a2      -> /compyfs/inputdata/atm/cam/chem/trop_mozart_aero/emis/DECK_ne30/cmip6_mam4_so4_a2_elev_1850-2014_c180205.nc'
ext_frc_type           = 'CYCLICAL'
ext_frc_cycle_yr       = 2010
srf_emis_specifier             = 'DMS       -> /compyfs/inputdata/atm/cam/chem/trop_mozart_aero/emis/DMSflux.2010.1deg_latlon_conserv.POPmonthlyClimFromACES4BGC_c20190220.nc',
        'SO2       -> /compyfs/inputdata/atm/cam/chem/trop_mozart_aero/emis/DECK_ne30/cmip6_mam4_so2_surf_1850-2014_c180205.nc',
        'bc_a4     -> /compyfs/inputdata/atm/cam/chem/trop_mozart_aero/emis/DECK_ne30/cmip6_mam4_bc_a4_surf_1850-2014_c180205.nc',
        'num_a1    -> /compyfs/inputdata/atm/cam/chem/trop_mozart_aero/emis/DECK_ne30/cmip6_mam4_num_a1_surf_1850-2014_c180205.nc',
        'num_a2    -> /compyfs/inputdata/atm/cam/chem/trop_mozart_aero/emis/DECK_ne30/cmip6_mam4_num_a2_surf_1850-2014_c180205.nc',
        'num_a4    -> /compyfs/inputdata/atm/cam/chem/trop_mozart_aero/emis/DECK_ne30/cmip6_mam4_num_a4_surf_1850-2014_c180205.nc',
        'pom_a4    -> /compyfs/inputdata/atm/cam/chem/trop_mozart_aero/emis/DECK_ne30/cmip6_mam4_pom_a4_surf_1850-2014_c180205.nc',
        'so4_a1    -> /compyfs/inputdata/atm/cam/chem/trop_mozart_aero/emis/DECK_ne30/cmip6_mam4_so4_a1_surf_1850-2014_c180205.nc',
        'so4_a2    -> /compyfs/inputdata/atm/cam/chem/trop_mozart_aero/emis/DECK_ne30/cmip6_mam4_so4_a2_surf_1850-2014_c180205.nc'
srf_emis_type          = 'CYCLICAL'
srf_emis_cycle_yr      = 2010
!
cflx_cpl_opt = 1
!
!...................
! conditional diag
!...................
metric_name = 'ALL',
!
qoi_chkpt = 'CFLX1','AERDRYRM','PBCINI','CFLX2',
            'CFLX3_01','CLDMAC01','CLDMIC01',
            'CFLX3_02','CLDMAC02','CLDMIC02',
            'CFLX3_03','CLDMAC03','CLDMIC03',
            'CFLX3_04','CLDMAC04','CLDMIC04',
            'CFLX3_05','CLDMAC05','CLDMIC05',
            'CFLX3_06','CLDMAC06','CLDMIC06',
            'AERWETRM',
!
qoi_name = 'dst_a1','dst_a1','dst_a3','dst_a3'
qoi_nver =  72,      72,      72,      72
qoi_x_dp =   0,       2,       0,       2
!
l_output_state = .true.
l_output_incrm = .true.
!
hist_tape_with_all_output = 1,2
!
/
EOF

cat <<EOB >! user_nl_clm
&clm_inparm
 finidat = '/compyfs/zhan524/myinput/init/compy_FC5AV1C-04P2_ne30_E3SM_20190418_ctrl_1600p_EXP01.clm2.r.0000-10-01-00000.nc'
 check_finidat_fsurdat_consistency = .false.
/
EOB

#####################################################################
# Job submission
#####################################################################

#./case.submit

endif

