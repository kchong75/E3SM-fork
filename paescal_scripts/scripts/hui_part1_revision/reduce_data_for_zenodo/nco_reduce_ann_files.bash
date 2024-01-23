###############################################################################################
# This script extracts variables from annually averaged history output.
#  - Input: history files in climo directories.
#  - Output: same file names but in different directories.
#
# Script history: Hui Wan (PNNL, 2023-12)
###############################################################################################
for cflx in 1 2; do
for nlev in 72 30; do

  casename="cflx_"$cflx"_L"$nlev"_F20TRC5-CMIP6"
  srcdir="/pscratch/sd/h/huiwan/cflx/v1_papers/v1_runs_original/"$casename"/climo/"
  dstdir="/pscratch/sd/h/huiwan/cflx/v1_papers/v1_runs_reduced/"$casename"/climo/"
  mkdir -p $dstdir
  
  yrS=2000
  yrE=2009
  dataFileTime="ANN"

  #=============================
  # multi-year mean annual mean
  #=============================
  yrRng=$yrS"-"$yrE
  filein=$casename".cam.h0."$yrRng"."$dataFileTime".nc"
  fileout=$filein

  echo $casename $yrRng

  #--------------------------------------------------------
  # part 1, fig03a,b: only need vertical integrals ("_v")
  #--------------------------------------------------------
  varlist="lon,lat,area"
  varlist=$varlist",cnd01_dst_a1_v_CFLX"$cflx"_inc,cnd01_dst_a3_v_CFLX"$cflx"_inc"
  varlist=$varlist",cnd01_dst_a1_v_AERWETRM,cnd01_dst_a3_v_AERWETRM"

  #-----------------------------------
  # part 2, fig02: need 3D tendencies
  #-----------------------------------
  varlist=$varlist",cnd01_dst_a1_CFLX"$cflx"_inc,cnd01_dst_a3_CFLX"$cflx"_inc"
  varlist=$varlist",cnd01_dst_a1_STCLD_inc,cnd01_dst_a3_STCLD_inc"
  
  ncks -v $varlist $srcdir/$filein $dstdir/$fileout

  echo
  #=============================
  # Yearly files
  #=============================
  year=$yrS
  while [ $year -le $yrE ]; do

    echo $casename $year

    filein=$casename".cam.h0."$year"."$dataFileTime".nc"
    fileout=$filein

    #---------------------------------------------
    # lon, lat, area, and emissions
    #---------------------------------------------
    varlist="lon,lat,area,cnd01_dst_a1_v_CFLX"$cflx"_inc,cnd01_dst_a3_v_CFLX"$cflx"_inc"

    #--------------------------------------
    # part 1, fig03, 07, 08: 3D tendencies
    #--------------------------------------
    for chkpt in "CFLX"$cflx "AERDRYRM" "PBCINI" "STCLD" "AERWETRM"; do
      varlist=$varlist",cnd01_dst_a1_"$chkpt"_inc"
      varlist=$varlist",cnd01_dst_a3_"$chkpt"_inc"
    done

    #--------------------------------------------------------------------------------------------
    # part 1, fig09: need burden and vertically integrated dry removal (in addition to emission)
    #--------------------------------------------------------------------------------------------
    varlist=$varlist",cnd01_dst_a1_v_AERWETRM,cnd01_dst_a3_v_AERWETRM"
    varlist=$varlist",cnd01_dst_a1_v_AERDRYRM_inc,cnd01_dst_a3_v_AERDRYRM_inc"

    #-----------------------------------------------------------------------------------------
    # part 1, table 1, vertical integrals: in addition to emissions, dry removal, burden,
    # also need STCLD (activation in this case), wet removal, and transport rates
    #-----------------------------------------------------------------------------------------
    varlist=$varlist",cnd01_dst_a1_v_STCLD_inc,cnd01_dst_a3_v_STCLD_inc"
    varlist=$varlist",cnd01_dst_a1_v_AERWETRM_inc,cnd01_dst_a3_v_AERWETRM_inc"
    varlist=$varlist",cnd01_dst_a1_v_PBCINI_inc,cnd01_dst_a3_v_PBCINI_inc"

    ncks -v $varlist $srcdir/$filein $dstdir/$fileout
    year=$(( $year + 1 ))
  done
  echo

done
done
