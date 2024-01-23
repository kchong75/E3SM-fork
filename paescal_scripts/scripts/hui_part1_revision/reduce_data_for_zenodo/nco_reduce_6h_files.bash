###############################################################################################
# This script copies 6-hourly output files but excludes unused variables to reduce file sizes.
#  - Input: history files in run directories.
#  - Output: same file names but in different directories.
#
# Script history: Hui Wan (PNNL, 2023-12)
###############################################################################################
for nlev in 72 30; do

  casename="cflx_1_L"$nlev"_F20TRC5-CMIP6_6h"
  srcdir="/pscratch/sd/h/huiwan/cflx/v1_papers/v1_runs_original/"$casename"/run/"
  dstdir="/pscratch/sd/h/huiwan/cflx/v1_papers/v1_runs_reduced/"$casename"/run/"
  mkdir -p $dstdir

  echo
  echo $casename

  #-----------------------------------------------------------------
  # We have multiple files to process for each casename (simulation)
  #-----------------------------------------------------------------
  for instOutputFileSuffix in \
  ".cam.h1.2000-01-29-00000.nc"  \
  ".cam.h1.2000-02-28-00000.nc"  \
  ".cam.h1.2000-03-30-00000.nc"
  do

    filein=$casename$instOutputFileSuffix
    fileout=$filein

    echo " Processing file "$filein

    # varlist is the list of variables to REMOVE from files

    varlist="cnd01_ALL_0e_to_360e_5n_to_55n,cnd01_ALL_flag_0e_to_360e_5n_to_55n"

    for mode in "a1" "a3"; do
        varlist=$varlist",cnd01_dst_"$mode"_CFLX2_0e_to_360e_5n_to_55n"
        varlist=$varlist",cnd01_dst_"$mode"_CFLX2_inc_0e_to_360e_5n_to_55n"
        varlist=$varlist",cnd01_dst_"$mode"_CHEM_0e_to_360e_5n_to_55n"
        varlist=$varlist",cnd01_dst_"$mode"_CHEM_inc_0e_to_360e_5n_to_55n"

        varlist=$varlist",cnd01_dst_"$mode"_v_AERDRYRM_0e_to_360e_5n_to_55n"
        varlist=$varlist",cnd01_dst_"$mode"_v_AERDRYRM_inc_0e_to_360e_5n_to_55n"
        varlist=$varlist",cnd01_dst_"$mode"_v_AERWETRM_0e_to_360e_5n_to_55n"
        varlist=$varlist",cnd01_dst_"$mode"_v_AERWETRM_inc_0e_to_360e_5n_to_55n"
        varlist=$varlist",cnd01_dst_"$mode"_v_CFLX1_0e_to_360e_5n_to_55n"
        varlist=$varlist",cnd01_dst_"$mode"_v_CFLX2_0e_to_360e_5n_to_55n"
        varlist=$varlist",cnd01_dst_"$mode"_v_CFLX2_inc_0e_to_360e_5n_to_55n"
        varlist=$varlist",cnd01_dst_"$mode"_v_CHEM_0e_to_360e_5n_to_55n"
        varlist=$varlist",cnd01_dst_"$mode"_v_CHEM_inc_0e_to_360e_5n_to_55n"
        varlist=$varlist",cnd01_dst_"$mode"_v_PBCINI_0e_to_360e_5n_to_55n"
        varlist=$varlist",cnd01_dst_"$mode"_v_PBCINI_inc_0e_to_360e_5n_to_55n"
        varlist=$varlist",cnd01_dst_"$mode"_v_STCLD_0e_to_360e_5n_to_55n"
        varlist=$varlist",cnd01_dst_"$mode"_v_STCLD_inc_0e_to_360e_5n_to_55n"
    done

    ncks -C -O -x -v $varlist $srcdir/$filein $dstdir/$fileout

  done #instOutputFileSuffix
done #nlev
