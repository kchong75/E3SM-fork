here=`pwd`
cd /pscratch/sd/h/huiwan/cflx/v1_papers/v1_runs_reduced/

for casename in \
    "cflx_1_L72_F20TRC5-CMIP6_6h" \
    "cflx_1_L30_F20TRC5-CMIP6_6h" \
do

  tar cvzf - ${casename}/ | split --bytes=3GB - ${casename}.tar.gz.

  #---------- mv the split tarballs to another dir and untar ------
  # mv ${casename}.tar.gz.* /pscratch/sd/h/huiwan/cflx/v1_papers/tarsplit/
  # cd /pscratch/sd/h/huiwan/cflx/v1_papers/tarsplit/
  # cat ${casename}.tar.gz.* | tar xzvf -
  #--------------------------------------

done

cd $here
