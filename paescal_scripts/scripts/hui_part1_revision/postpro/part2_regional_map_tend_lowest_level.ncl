  load "./utils/contour_line_utils.ncl"

  cvd_friendly  = True
  cvd_friendly  = False

  itime   = 0
  year    = "2000-2009"
  dtime    = 1800.
  nlevEAM  = 72

  expname = (/"cflx_1_L72_F20TRC5-CMIP6",\
              "cflx_2_L72_F20TRC5-CMIP6"/)
  expshortname = (/"v1_ori","v1_rev"/)+"_L"+nlevEAM
  expshortname = (/"original","revised"/)
  nexp = dimsizes(expname)

  climoFilePath = "/pscratch/sd/h/huiwan/cflx/v1_papers/v1_runs/"+expname+"/climo/"
  FileName = expname+".cam.h0."+year+".ANN.nc"

  plotPath = "/global/cfs/projectdirs/m4359/www/huiwan/cflx/2023_rerun_for_gmd/"
  plotFile = "dust_A_C."+year
  plotFmt  = "pdf"

  label = (/"Emissions","Turb-mix+act/resus"/)
  chkpt = (/(/"CFLX1","STCLD"/),\
            (/"CFLX2","STCLD"/) /)
  QoI_name = "cnd01_"+(/"dst_a1","dst_a3"/)
  inc_fac  = 1./dtime

  v_suff   = ""                         ; not vertical integrals
  inc_suff = "_inc"                     ; "_inc" means increment; dividing by dtime gives tendency
  unit     = "kg kg~S~-1~N~ s~S~-1~N~"

  ;-----------------
  ; What to plot
  ;-----------------
  term_coeff_for_panels = (/  \
                             (/ 1.,  0./), \    ;  A+0
                             (/ 0.,  1./), \    ;  0+C
                             (/ 1., -1./), \    ;  A-C
                             (/-1., -1./)  /)   ; -A-C

  panel_title  = (/"A","C", "A-C","-A-C"/)
  panel_prefix = (/ (/"","", "(a)","(b)"/), \
                    (/"","", "(c)","(d)"/) /)
  ;-----------------
  ; Colormaps
  ;-----------------

  if (cvd_friendly) then
     plotFile = plotFile +"_cvd"
     cmap = read_colormap_file("bam")
     diff_colormap  = cmap(::-1,:)     ; reverse color map
     delete(cmap)
    ;cnLevels = (/-1e-8,-3e-9, -1e-9,-3e-10,-1e-10,1e-10,3e-10,1e-9,3e-9,1e-8/) 
     cnLevels = (/ -5e-9, -2e-9, -1e-9,-5e-10,-2e-10,-1e-10,1e-10,2e-10,5e-10,1e-9,2e-9,5e-9/) 
  else
     plotFile = plotFile
     diff_colormap = "BlueWhiteOrangeRed"
     cnLevels = (/ -5e-9, -2e-9, -1e-9,-5e-10,-2e-10,-1e-10,1e-10,2e-10,5e-10,1e-9,2e-9,5e-9/) 
  end if
  ;-----------------


  nchkpt = dimsizes(label)
  nqoi   = dimsizes(QoI_name)

  npanel = dimsizes(panel_title)

  ;-----------------------------
  ; for plotting
  ;-----------------------------
  wks  = gsn_open_wks(plotFmt,plotPath+plotFile)
  system("chmod a+r "+plotPath+plotFile+"."+plotFmt)

  res = True
  res@gsnFrame = False
  res@gsnDraw  = False

 ;res@mpProjection         = "Robinson"
 ;res@mpDataBaseVersion    = "LowRes"       ; better map outlines
 ;res@mpGridAndLimbOn      = True           ; turn on limb and grid lines
 ;res@mpGridLineColor      = -1             ; ...but don't draw grid lines...
 ;res@mpPerimOn            = False          ; turn off map perimeter
 ;res@mpFillDrawOrder      = "PostDraw"     ; draw map fill last
 ;res@mpGeophysicalLineColor = "gray"

  res@mpLimitMode = "LatLon"
  res@mpMinLatF =  5
  res@mpMaxLatF = 40           ; choose subregion           
  res@mpMinLonF = -20
  res@mpMaxLonF =  35

  res@mpFillOn    = False                       ; turn map fill off
  res@mpOutlineOn = True                        ; turn the map outline on
  res@cnFillOn    = True                        ; turn on color fill
  res@cnLineLabelsOn = False                    ; turn the line labels off

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

  res@cnFillPalette  = diff_colormap
  ;-----------------------------
  ; open data file
  ;-----------------------------
  File = addfiles(climoFilePath+FileName,"r")

  area = tofloat(File[0]->area)
  lat  = File[0]->lat
  lon  = File[0]->lon
; lon  = where(lon.le.180., lon, lon - 360.)

  res@sfXArray = lon
  res@sfYArray = lat

  res@cnLevelSelectionMode = "ExplicitLevels"   ; set explicit contour levels
  ;---------------------


plot = new( (/nexp,npanel/),"graphic" )

do iexp = 0,nexp-1

  process_rates = new( (/nchkpt,dimsizes(lat)/), "float" )
  process_rates = 0. 

  do ichkpt = 0,nchkpt-1 ; loop over different processes
     do iqoi = 0,nqoi-1 ; for each process, sum over all modes
        varname = QoI_name(iqoi)+v_suff+"_"+ chkpt(iexp,ichkpt)+inc_suff
        process_rates(ichkpt,:) = process_rates(ichkpt,:) + File[iexp]->$varname$(itime,nlevEAM-1,:)*inc_fac
     end do
  end do
  print("Done reading  data for "+expname(iexp))

  do ipanel=0,npanel-1

     var2d = new( dimsizes(lat), "float")
     var2d = 0.
     do ichkpt=0,nchkpt-1
        var2d = var2d + process_rates(ichkpt,:)*term_coeff_for_panels(ipanel,ichkpt)
     end do

     res@tiMainString = panel_prefix(iexp,ipanel)+" "+panel_title(ipanel) +" in EAMv1 with "+expshortname(iexp)+" coupling"
     res@cnLevels = cnLevels 

     plot(iexp,ipanel) = gsn_csm_contour_map(wks,var2d,res)    

     ; Use dashed lines for contours of negative values, for being CVD-friendly.
     ; (In this particular script, because we are using the same contour levels and colormap for all panels,
     ; it is unnecessary to assign res@cnLevels, cnDashPatterns, and cnThicknesses and delete them repeatedly.
     ; I'm keeping the corresponding lines here just in case they are needed in future examples.)
    ;if (cvd_friendly) then

    ;   res@cnLinesOn   = True 
    ;   cnDashPatterns = toint(where(res@cnLevels.lt.0, 2,  0  ))
    ;   cnThicknesses  =       where(res@cnLevels.lt.0, 1., 0.1)
    ;   set_contour_line_patterns   (plot(iexp,ipanel),res@cnLevels,cnDashPatterns)
    ;   set_contour_line_thicknesses(plot(iexp,ipanel),res@cnLevels,cnThicknesses)

    ;   delete([/cnDashPatterns,cnThicknesses/])
    ;   delete(res@cnLevels)
    ;else
        res@cnLinesOn   = False 
    ;end if

  end do
  print("Done plotting data for "+expname(iexp))

  resP = True
  resP@gsnPanelXWhiteSpacePercent = 5 
  resP@gsnPanelYWhiteSpacePercent = 5 
  resP@gsnPanelMainFontHeightF = res@gsnStringFontHeightF*0.3
 ;gsn_panel(wks,plot(iexp,:),(/1,npanel/),resP)

end do

  gsn_panel(wks,(/plot(0,2),plot(0,3),plot(1,2),plot(1,3)/),(/2,2/),resP)

  print("plot file is "+(plotPath+plotFile+"."+plotFmt))
