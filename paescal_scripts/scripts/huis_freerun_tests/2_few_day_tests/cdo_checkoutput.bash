histfile="/compyfs/wanh895/scidac_scratch/2023_cflx/2023code/opt40_F20TRC5-CMIP6/run/case.cam.h0.2009-10-03-00000.nc"
histfile="/compyfs/wanh895/scidac_scratch/2023_cflx/2023code/opt41_F20TRC5-CMIP6/run/case.cam.h0.2009-10-03-00000.nc"

cdo infov -selname,dst_a1DDF,dst_a3DDF $histfile 
cdo infov -selname,dst_a1TBF,dst_a3TBF $histfile
cdo infov -selname,dst_a1GVF,dst_a3GVF $histfile

#cdo infov -selname,num_a1_TBV,mss_a1_TBV $histfile
#cdo infov -selname,num_a3_TBV,mss_a3_TBV $histfile
#cdo infov -selname,num_a1_GVV,mss_a3_GVV $histfile
#cdo infov -selname,num_a3_GVV,mss_a3_GVV $histfile

