  load "./utils/sub_lonlat_utils.ncl" ; utility functions

  l_debug = True
  l_debug = False

  ;------------
  ; Time range
  ;------------
  yrS   = 2000
  yrE   = 2009
  nyr   = yrE - yrS + 1

  yrRng  = yrS+"-"+yrE
  simlen = nyr+" years"

  dataFileTime = "ANN"

  ;-----------
  ; casenames
  ;------------
  casename = (/"cflx_1_L72_F20TRC5-CMIP6", \
               "cflx_2_L72_F20TRC5-CMIP6"  /)

  ncase = dimsizes(casename)

  dataFilePath = "/pscratch/sd/h/huiwan/cflx/v1_papers/v1_runs/"+casename+"/climo/"

  itime = 0      ; each file contains only one time slice
  dtime = 1800.  ; model stepsize, needed for converting increments to tendencies

  ;------------
  ; Plot file

  plotPath = "/global/cfs/projectdirs/m4359/www/huiwan/cflx/2023_rerun_for_gmd/"
  plotFileTimeTag = yrRng+"_"+dataFileTime

  ;------------------------------------------------------------------------
  ; Aerosol species and the corresponding mixing ratios in different modes
  ;------------------------------------------------------------------------
  species  = "dst"
  QoI_name = "cnd01_"+(/"dst_a3","dst_a1"/)
  nqoi     = dimsizes(QoI_name)

  ;--------------------------------------------------
  ; Checkpoints
  ;--------------------------------------------------
  chkpt    = (/(/"CFLX1","AERDRYRM","PBCINI", "STCLD","AERWETRM"/),\
               (/"CFLX2","AERDRYRM","PBCINI", "STCLD","AERWETRM"/)/)
  label    = (/"Emissions","Dry removal","Resolved transport","Turb-mix+act/resus","Wet removal"/)
  nchkpt   = dimsizes(label)

  ichkpt_emis = ind(label.eq."Emissions")
  ichkpt_wtrm = ind(chkpt(0,:).eq."AERWETRM")

  print("")
  if (ichkpt_emis.lt.0) then
     print("Error: did not find checkpoint corresponding to emissions - needed for conditional average")
     exit
  else
     print("casename = "+casename+": chkpt corresponding to emissions = "+(ichkpt_emis)+": "+chkpt(:,ichkpt_emis))
     print("")
  end if

  ;------------------------------------------------------------------------------
  ; Regions in which the mean vertical profiles will be calculated and plotted
  ;------------------------------------------------------------------------------
  region_label = (/"Source-region","Vicinity-region","Remote-region"/)

  latmin    = (/  5,  5, 10/)
  latmax    = (/ 55, 55, 40/)
  lonmin    = (/-25,-25,160/)
  lonmax    = (/130,130,230/)
  emis_type = (/  1,  0,  0/)     ; 1 = has emission; 0 = no emission

  nregion = dimsizes(emis_type)


  ;=======================================================================================================
  FontHeightF = 0.02

  ; The following settings are used both for the global maps of emission/burden and for
  ; sanity checks before calculting regional averages, and hence are placed here in the main/driver script.

  cnres = True
  cnres@gsnFrame = False
  cnres@gsnDraw  = False

  cnres@mpProjection         = "Robinson"
  cnres@mpDataBaseVersion    = "LowRes"       ; better map outlines
  cnres@mpGridAndLimbOn      = True           ; turn on limb and grid lines
  cnres@mpGridLineColor      = -1             ; ...but don't draw grid lines...
  cnres@mpPerimOn            = False          ; turn off map perimeter
  cnres@mpFillDrawOrder      = "PostDraw"     ; draw map fill last
  cnres@mpGeophysicalLineColor = "gray70"     ; light gray. Larger numbers are lighter.
  cnres@mpCenterLonF           = 155. + 180.

  cnres@mpFillOn    = False                       ; turn map fill off
  cnres@mpOutlineOn = True                        ; turn the map outline on
  cnres@cnFillOn    = True                        ; turn on color fill
  cnres@cnLinesOn   = False                       ; turn off the contour lines
  cnres@cnLineLabelsOn = False                    ; turn the line labels off

  cnres@tmXTOn = False
  cnres@tmYROn = False

  cnres@gsnLeftString = ""
  cnres@gsnCenterString = ""
  cnres@gsnRightString = ""

  cnres@tiMainFontHeightF    = FontHeightF*1.1
  cnres@gsnStringFontHeightF = FontHeightF

 ;cnres@lbOrientation = "horizontal"
  cnres@lbOrientation = "vertical"
  cnres@pmLabelBarHeightF  = 0.41
  cnres@lbLabelFontHeightF = cnres@gsnStringFontHeightF

  ;=======================================================================================================

  load "./fig03_global_maps.ncl"  ; plot global maps of emissions and burden. Add boxes to the emission map

  load "./fig03_08_regional_mean_tend_profiles.ncl" ; calculate and plot regional mean profiles.

  ;=======================================================================================================
