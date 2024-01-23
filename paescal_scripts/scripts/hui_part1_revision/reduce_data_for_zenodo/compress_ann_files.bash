  srcdir="/pscratch/sd/h/huiwan/cflx/v1_papers/v1_runs_reduced/"
  cd $srcdir

  for cflx in 1 2; do
  for nlev in 72 30; do

      casename="cflx_"$cflx"_L"$nlev"_F20TRC5-CMIP6"
      tar -czvf $casename".tar.gz" $casename

  done
  done
