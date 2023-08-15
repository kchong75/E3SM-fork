  ;-----------------------------------------------------------
  nlevEAM = 30
  nlevEAM = 72

  if (nlevEAM.eq.72) then
     casename = (/"compy_F20TRC5-CMIP6_ne30_ndg_cflx_cpl_opt1",\
                  "compy_F20TRC5-CMIP6_ne30_ndg_cflx_cpl_opt2" /)
  else
     casename = (/"compy_F20TRC5-CMIP6_ne30_ndg_cflx_cpl_opt1_L"+nlevEAM,\
                  "compy_F20TRC5-CMIP6_ne30_ndg_cflx_cpl_opt2_L"+nlevEAM /)
  end if

  climFilePath = "/compyfs/wanh895/scidac_scratch/2023_cflx/2021code_kai/"+casename+"/climo/"
  histFilePath = "/qfs/people/zhan524/cshr/e3sm_scratch/cflx/simulation/"+casename+"/run/"

  plotFileSuffix = "."+casename

  ;------------------------------------------------------
  ; Aerosol species and the corresponding mixing ratios
  ;------------------------------------------------------
  species  = "dst"
  QoI_name = "cnd01_"+(/"dst_a3","dst_a1"/)
  nqoi     = dimsizes(QoI_name)

  dtime     = 1800.
  year      = 2010
  remap_res = "180x360"
  ;--------------------------------------------------
  ; Checkpoints
  ;--------------------------------------------------
  chkpt    = (/(/"CFLX1","AERDRYRM","PBCINI", "STCLD","AERWETRM"/),\
               (/"CFLX2","AERDRYRM","PBCINI", "STCLD","AERWETRM"/)/)
  label    = (/"Emissions","Dry removal","Resolved transport","Turb-mix+act/resus","Wet removal"/)
  nchkpt   = dimsizes(label)

  ichkpt_emis = ind(label.eq."Emissions")

  print("")
  if (ichkpt_emis.lt.0) then
     print("Error: did not find checkpoint corresponding to emissions - needed for conditional average")
     exit
  else
     print("chkpt corresponding to emissions = "+(ichkpt_emis+1))
  end if

  ;=======================================================================================================
  ; Regions in which the mean vertical profiles will be calculated and plotted
  ;=======================================================================================================
  latmin = (/  5,  5, 10/)
  latmax = (/ 55, 55, 40/)
  lonmin = (/-25,-25,160/)
  lonmax = (/130,130,230/)
  emis_type = (/1,0,0/)     ; 1 = with emission; 0 = no emission

  nregion = dimsizes(emis_type)

  region_label = (/"Source-region","Vicinity-region","Remote-region"/)
  region_panel_lb = (/"(a)","(b)","(c)"/)
  region_panel_lb = "" 


  tmp_data_dir = "./data"+"_"+casename+"/"

  do icase = 0,dimsizes(casename)-1
     system("mkdir -p "+tmp_data_dir(icase))
  end do

  plotFilePath = "./"
