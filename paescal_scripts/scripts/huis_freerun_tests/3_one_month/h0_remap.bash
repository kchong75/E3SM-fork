#! /bin/bash 

module load cdo nco

hist_root="/compyfs/wanh895/scidac_scratch/2023_cflx/2023code/1mon/refac/"
clim_root="/compyfs/wanh895/scidac_scratch/2023_cflx/2023code/1mon/refac/"
year=2009
month=10

hX="h0"

res_latlon="180x360"


#opt41_F20TRC5-CMIP6 \
#opt42_F20TRC5-CMIP6 \
#opt43_F20TRC5-CMIP6 \
#opt44_F20TRC5-CMIP6 
for casename in \
opt1_F20TRC5-CMIP6 \
opt2_F20TRC5-CMIP6 \
opt3_F20TRC5-CMIP6 \
opt4_F20TRC5-CMIP6 
do

  cd ${hist_root}/${casename}/run
  echo
  echo "Processing "${hX}" output in "`pwd`
  echo
  mkdir -p ${clim_root}/${casename}/remap

  ls ${casename}.cam.${hX}.${year}-??.nc
  echo
  echo

  if [[ $hX == "h0" ]]; then
         ncremap -m "/qfs/people/zender/data/maps/map_ne30np4_to_cmip6_"${res_latlon}"_aave.20181001.nc" \
             ${clim_root}/${casename}/run/${casename}.cam.${hX}.${year}-${month}.nc \
             ${clim_root}/${casename}/remap/${casename}.cam.${hX}.${year}-${month}.${res_latlon}.nc
  fi

done
