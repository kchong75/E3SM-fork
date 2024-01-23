  load "./utils/contour_line_utils.ncl"
 
  l_debug = True
  l_debug = False

  FillValue = -999.

  ;------------
  ; Time range
  ;------------
  yrS   = 2000
  yrE   = 2009
  nyr   = yrE - yrS + 1

  yrRng  = yrS+"-"+yrE
  simlen = nyr+" years"

  dataFileTime = "ANN"

  itime   = 0

  ;-----------------------------
  ; casenames etc.
  ;-----------------------------
  dtime    = 1800.
  nlevEAM  = 72

  casename = (/"cflx_1_L72_F20TRC5-CMIP6", \
               "cflx_2_L72_F20TRC5-CMIP6"  /)

  caseShortName = (/"original","revised"/)
  ncase = dimsizes(casename)

  dataFilePath = "/pscratch/sd/h/huiwan/cflx/v1_papers/v1_runs/"+casename+"/climo/"

  label = (/"Emissions","Dry removal"/)
  chkpt = (/(/"CFLX1","AERDRYRM"/),\
            (/"CFLX2","AERDRYRM"/) /)

  species  = "dst"
  QoI_name = "cnd01_"+(/"dst_a1","dst_a3"/)

 ;cnLevels = (/ -3e-10, -1e-10, -3e-11,-1e-11,-3e-12, 3e-12,1e-11,3e-11, 1e-10,  3e-10/) ; more levels to show more detailed structures 
  cnLevels = (/ -1e-9,-3e-10, -1e-10, -3e-11,0.,3e-11, 1e-10,  3e-10,1e-9/)              ; fewer levels to accommodate CVD-friendly colormaps

  nchkpt = dimsizes(chkpt(0,:))
  nqoi   = dimsizes(QoI_name)

      unit = "kg kg~S~-1~N~ s~S~-1~N~"
    v_suff = ""                             ; 3D field from CondiDiag, no vertical integral
  inc_suff = "_inc"                         ; increments from CondiDiag 
  inc_fac  = 1./dtime                       ; factor for converting increments to tendencies

  ;-----------------------------
  ; for plotting
  ;-----------------------------
  panel_title_prefix = (/ (/"(a)","   ","   "/), \
                          (/"(b)","(c)","(d)"/)  /)

  plot  = new( (/nchkpt,3/),"graphic" )

  plotFileTimeTag = yrRng+"_"+dataFileTime
  plotPath = "/global/cfs/projectdirs/m4359/www/huiwan/cflx/2023_rerun_for_gmd/"
  plotFile = plotPath+species+"_emis_dryrm_tend_lowest_layer_cmpr_"+plotFileTimeTag
  plotFmt  = "pdf"

  wks      = gsn_open_wks(plotFmt,plotFile)
  system("chmod a+r "+plotFile+"."+plotFmt)
  gsn_define_colormap(wks,"bam_with_white")
  gsn_reverse_colormap(wks)

  print("Plot file is "+plotFile+"."+plotFmt)

  res = True
  res@gsnFrame = False
  res@gsnDraw  = False

  res@mpLimitMode = "LatLon"
  res@mpMinLatF =  5
  res@mpMaxLatF = 40           ; choose subregion           
  res@mpMinLonF = -20
  res@mpMaxLonF =  35

  res@mpFillOn    = False                   ; turn map fill off
  res@mpOutlineOn = True                    ; turn the map outline on
  res@mpGeophysicalLineColor = "gray60"     ; light gray. Larger numbers are lighter.


  res@cnFillOn    = True
  res@cnLinesOn   = True
  res@cnLineLabelsOn = False 

  res@tmXTOn = False
  res@tmYROn = False

  FontHeightF = 0.03
  res@tiMainFontHeightF    = FontHeightF*1.
  res@gsnStringFontHeightF = FontHeightF
  res@tiXAxisFontHeightF   = FontHeightF*0.8
  res@tiYAxisFontHeightF   = FontHeightF*0.8
  res@tmXBLabelFontHeightF = FontHeightF*0.6
  res@tmYLLabelFontHeightF = FontHeightF*0.6

  res@lbOrientation = "vertical"
  res@lbLabelFontHeightF = res@gsnStringFontHeightF * 0.8
  res@pmLabelBarHeightF   = 0.5
 ;res@pmLabelBarWidthF    = 0.1
 ;res@pmLabelBarOrthogonalPosF = 0.15
 ;res@pmLabelBarParallelPosF = 0.6

  ;---------------------------------------
  ; Open a data file and get grid info
  ;--------------------------------------
  dataFileName = casename(0)+".cam.h0."+yrS+"."+dataFileTime+".nc"
  dataFile = addfile(dataFilePath(0)+dataFileName,"r")

  lat  = dataFile->lat
  lon  = dataFile->lon
  ncol = dimsizes(lon)

  res@sfXArray = lon
  res@sfYArray = lat

  ;---------------------
  var2d = new((/ncase,ncol,nyr/),"float")

  do ichkpt = 0,nchkpt-1

    print("")
    print(" Process "+ichkpt)
    print("")
    ;-----------------------------------------------------------
    ; Read in model output and sum over different aerosol modes
    ;-----------------------------------------------------------
    do icase = 0,ncase-1
    do iyr = 0,nyr-1

       ; For each case and each year we need to read data from a different file

       dataFileName = casename(icase)+".cam.h0."+(yrS+iyr)+"."+dataFileTime+".nc"
       dataFile = addfile(dataFilePath(icase)+dataFileName,"r")

       ; Initialize, then calculate sum.

       var2d(icase,:,iyr) = 0.

       do iqoi = 0,nqoi-1
          varname = QoI_name(iqoi)+v_suff+"_"+ chkpt(icase,ichkpt)+inc_suff
          var2d(icase,:,iyr) = var2d(icase,:,iyr) + dataFile->$varname$(itime,nlevEAM-1,:)*inc_fac
       end do

       print("Done reading data from "+dataFileName)
     end do
     print("")
     end do

     ;-------------------------------------------
     ; Calculate multi-year mean and std. dev.
     ;------------------------------------------
     var2d_mean = dim_avg   ( var2d )
     var2d_std  = dim_stddev( var2d ) 

     ;------------------------------------------
     ; For plotting
     ;------------------------------------------
     res@cnLevelSelectionMode = "ExplicitLevels"   ; set explicit contour levels
     res@cnLevels = cnLevels

     cnDashPatterns = toint(where(res@cnLevels.lt.0, 2,  0  ))
     cnThicknesses  =       where(res@cnLevels.lt.0, 1., 0.1)

     ;-------------
     ; control
     ;-------------
     ii = 0     ; first column in panel plot

     icase = 0  ; control case
     res@tiMainString = panel_title_prefix(ichkpt,ii)+" "+label(ichkpt) +" in layer "+nlevEAM
     res@gsnCenterString = "("+caseShortName(icase)+", unit: "+unit+")"

     plot(ichkpt,ii) = gsn_csm_contour_map(wks,var2d_mean(icase,:),res)    

     set_contour_line_patterns   (plot(ichkpt,ii),res@cnLevels,cnDashPatterns)
     set_contour_line_thicknesses(plot(ichkpt,ii),res@cnLevels,cnThicknesses)

     ;----------------------------------------------
     ; Perform a Kolmogorov-Smirnov two-sample test 
     ;----------------------------------------------
     probKS = kolsm2_n(var2d(1,:,:), var2d(0,:,:), 1)  ; KS test is applied to the rightmost dimension
     signif = 0.05

     ;-------------
     ; diff.
     ;-------------
     ii = 1     ; second column in panel plot

     res@tiMainString    = panel_title_prefix(ichkpt,ii)+" "+label(ichkpt) +" in layer "+nlevEAM+ ", difference"
     res@gsnCenterString = "("+caseShortName(1)+" vs. "+caseShortName(0)+", unit: "+unit+")"

     diff = where( probKS.lt.signif, var2d_mean(1,:)-var2d_mean(0,:), FillValue )
     diff@_FillValue = FillValue

     plot(ichkpt,ii) = gsn_csm_contour_map(wks,diff,res)    
  
     set_contour_line_patterns   (plot(ichkpt,ii),res@cnLevels,cnDashPatterns)
     set_contour_line_thicknesses(plot(ichkpt,ii),res@cnLevels,cnThicknesses)

     ;------------------
     ; relative diff.
     ;------------------
     ii = 2    ; 3rd column in panel plot

     res@tiMainString    = panel_title_prefix(ichkpt,ii)+" "+label(ichkpt) +" in layer "+nlevEAM+ ", relative diff."
     res@gsnCenterString = "("+caseShortName(1)+" vs. "+caseShortName(0)+", unit: -)"

     delete(res@cnLevels)
     res@cnLevels = (/-1.5,-1.0,-0.50,0.,0.50,1.0,1.5/) 

     var_eps = 1e-36
     denom = where(abs(var2d_mean(0,:)).ge.var_eps, var2d_mean(0,:), FillValue)
     denom@_FillValue = FillValue 

     reldiff = where( probKS.lt.signif, (var2d_mean(1,:)-var2d_mean(0,:))/denom, FillValue )
     reldiff@_FillValue = FillValue

     plot(ichkpt,ii) = gsn_csm_contour_map(wks,reldiff,res)    

     delete([/cnDashPatterns,cnThicknesses/])
     cnDashPatterns = toint(where(res@cnLevels.lt.0, 2,  0  ))
     cnThicknesses  =       where(res@cnLevels.lt.0, 1., 0.1)

     set_contour_line_patterns   (plot(ichkpt,ii),res@cnLevels,cnDashPatterns)
     set_contour_line_thicknesses(plot(ichkpt,ii),res@cnLevels,cnThicknesses)

     delete(res@cnLevels)
     delete([/cnDashPatterns,cnThicknesses/])

  end do ; checkpoint loop

  resP = True
  resP@gsnMaximize = True
  resP@gsnPanelYWhiteSpacePercent = 8 
  resP@gsnPanelMainFontHeightF = res@gsnStringFontHeightF*0.3

  gsn_panel(wks,(/plot(0,0),plot(1,0),plot(1,1),plot(1,2)/),(/1,4/),resP)

  if (l_debug) then
     gsn_panel(wks,ndtooned(plot),(/nchkpt,3/),resP)
  end if

