###! /bin/bash 

process_single_year_ann=1
process_multi_year_mean=1

yearS=2000
yearE=2009

yrRng="200[0-9]"  # for cdo command that creates multi-year mean

hX="h0"

hist_root="/pscratch/sd/h/huiwan/cflx/v1_papers/v1_runs/"
clim_root="/pscratch/sd/h/huiwan/cflx/v1_papers/v1_runs/"

remap_wgtfile_path="/global/homes/z/zender/data/maps/"
res_latlon="180x360"

#=========================================
# Process one case (simulation) at a time
#=========================================
for casename in \
cflx_1_L30_F20TRC5-CMIP6 \
cflx_2_L30_F20TRC5-CMIP6 \
cflx_1_L72_F20TRC5-CMIP6 \
cflx_2_L72_F20TRC5-CMIP6
do

  echo
  echo "Processing "${hX}" output in "${clim_root}/${casename}
  echo

  mkdir -p ${clim_root}/${casename}/climo

  year=$yearS

  #----------------------------
  # Process one year at a time
  #----------------------------
  if [ $process_single_year_ann -ge 1 ]; then

     cd ${hist_root}/${casename}/run

     while [ $year -le $yearE ]
     do
        echo
        ls ${casename}.cam.${hX}.${year}-??.nc
        echo

        #------------------------------------------------------------------
        # Calculate and write out annual mean of single year on model grid
        #------------------------------------------------------------------
        cdo ensavg ${casename}.cam.${hX}.${year}-??.nc \
                   ${clim_root}/${casename}/climo/${casename}.cam.${hX}.${year}.ANN.nc
        echo
        echo Done calculating annual mean.
        echo

        #----------------------------------
        # Remap annual mean to latlon grid 
        #----------------------------------
        ncremap -m ${remap_wgtfile_path}"/map_ne30np4_to_cmip6_"${res_latlon}"_aave.20181001.nc" \
                ${clim_root}/${casename}/climo/${casename}.cam.${hX}.${year}.ANN.nc \
                ${clim_root}/${casename}/climo/${casename}.cam.${hX}.${year}.ANN.${res_latlon}.nc

        echo
        echo Done remapping annual mean to latlon grid.
        echo

        year=$(( $year + 1 ))
     done
  fi

  #--------------------------------------------------------
  # All individial years processed. Obtain multi-year mean
  #--------------------------------------------------------
  if [ $process_multi_year_mean -ge 1 ]; then

     cd ${clim_root}/${casename}/climo/

     cdo ensavg ${casename}.cam.${hX}.${yrRng}.ANN.nc \
                ${casename}.cam.${hX}.${yearS}-${yearE}.ANN.nc

     cdo ensavg ${casename}.cam.${hX}.${yrRng}.ANN.${res_latlon}.nc \
                ${casename}.cam.${hX}.${yearS}-${yearE}.ANN.${res_latlon}.nc

     echo
     echo Done processing case ${casename}
     echo
  fi

done # case loop
#=========================================
