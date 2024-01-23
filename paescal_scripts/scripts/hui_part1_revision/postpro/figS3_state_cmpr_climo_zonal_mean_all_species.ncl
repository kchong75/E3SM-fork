;---------------------------------------------------------------------------------------
; This script compares the zonal mean mass mixig ratios from two simulations.
; Different rows are different species.
; Column 1 shows the results from the second case (simulation);
; Column 2 shows the differences between the second and the first case;
; Column 3 shows the relative differences.
;---------------------------------------------------------------------------------------
; History:
;  - Initial version by Hui Wan, PNNL, spring 2022.
;  - Revision by Hui Wan, PNNL, 2023-11, to make the plots CVD-friendly.
;  - Revision by Hui Wan, PNNL, 2023-11, add option to plot cloudborne + interstitial
;    instead of interstitial aerosols only. 
;---------------------------------------------------------------------------------------
  load "./utils/contour_line_utils.ncl"

  cvd_friendly  = False
  cvd_friendly  = True

  l_debug = True
  l_debug = False

  FillValue = -999.           
  signif = 0.05      ; significance level used in statitistical tests

  ;------------
  ; Time range
  ;------------
  yrS   = 2000
  yrE   = 2009
  nyr   = yrE - yrS + 1

  yrRng  = yrS+"-"+yrE
  simlen = nyr+" years"

  dataFileTime = "ANN"

  ;-----------------------------
  ; casenames etc.
  ;-----------------------------

  casename = (/"cflx_1_L72_F20TRC5-CMIP6", \
               "cflx_2_L72_F20TRC5-CMIP6"  /)

 ;casename = (/"cflx_1_L30_F20TRC5-CMIP6", \
 ;             "cflx_2_L30_F20TRC5-CMIP6"  /)

  dataFilePath = "/pscratch/sd/h/huiwan/cflx/v1_papers/v1_runs/"+casename+"/climo/"

  species_category = 1 ; gases 
  species_category = 2 ; aerosols (interstitial only)
  species_category = 3 ; aerosols (interstitial + cloudborne)

  if (species_category.eq.1) then
     species_longname = (/"O3","H2O2","H2SO4","SO2","DMS","SOAG"/)
     species          = species_longname
     plotFileName     = "state_cmpr_all_gases_"+casename(1)+"-vs-"+casename(0)
  else
  if (species_category.eq.2 .or. species_category.eq.3) then
     species_longname = (/"Dust","Sea salt","MOA","BC","POA"/)
     species          = (/"dst", "ncl",     "mom","bc","pom"/)
     plotFileName     = "state_cmpr_all_aerosols_"+species_category+"_"+casename(1)+"-vs-"+casename(0)
  else
     print("Unknown species_category. Abort.")
  end if
  end if

  ;-------------
  ; Color maps
  ;-------------
   if (.not.cvd_friendly) then  ; use NCL colormaps

     state_colormap = "WhiteBlueGreenYellowRed"
      diff_colormap = "nrl_sirkes"

   else ; use CVD-friendly colormaps

     cmap = read_colormap_file("batlowW")
     state_colormap = cmap(::-1,:)
     delete(cmap)

     cmap = read_colormap_file("bam")
     diff_colormap  = cmap(::-1,:)     ; reverse color map
     delete(cmap)

   end if

  ;-----------------------------
  ; For plotting
  ;-----------------------------
  ncol = 3  ; (1) revised model, (2) revised - original, (3) relative diff.

  nrow = dimsizes(species)
  plot = new( (/nrow,ncol/),"graphic" )
   dum = new( (/nrow,ncol/),"graphic" )

  plotFileTimeTag = yrRng+"_"+dataFileTime

  plotPath = "/global/cfs/projectdirs/m4359/www/huiwan/cflx/2023_rerun_for_gmd/"
  plotFile = plotPath+plotFileName+"_"+plotFileTimeTag
  plotFmt  = "pdf"


  wks      = gsn_open_wks(plotFmt,plotFile)
  system("chmod a+r "+plotFile+"."+plotFmt)

  res = True
  res@gsnFrame = False
  res@gsnDraw  = False

  res@cnFillOn       = True
  res@cnLineLabelsOn = False

  res@tmXTOn = False
  res@tmYROn = False

  FontHeightF = 0.02 *0.35
  res@gsnStringFontHeightF = FontHeightF *0.95
  res@tiXAxisFontHeightF   = FontHeightF *0.95
  res@tiYAxisFontHeightF   = FontHeightF *0.95

  res@lbOrientation        = "vertical"
  res@lbLabelFontHeightF   = FontHeightF *0.8

  vpWthF_crs = 0.18
  vpHgtF_crs = vpWthF_crs *0.65
  res@vpWidthF  = vpWthF_crs 
  res@vpHeightF = vpHgtF_crs

  res@gsnLeftStringOrthogonalPosF = 0.035
  res@gsnRightStringOrthogonalPosF = 0.035

  ;-----------------------------------------------------
  ; For stippling that masks out insignificant results
  ;-----------------------------------------------------
  gsn_define_colormap(wks,"GMT_gray")

  stiplres                      = True             ; plot mods desired
  stiplres@gsnFrame             = False
  stiplres@gsnDraw              = False
  stiplres@tmXTOn               = False
  stiplres@tmYROn               = False
  stiplres@cnFillOn             = True             ; fill contour intervals
  stiplres@cnLinesOn            = False
  stiplres@cnLineLabelsOn       = False
  stiplres@cnInfoLabelOn        = False

  stiplres@cnLevelSelectionMode = "ManualLevels"   ; manually specify contour levels
  stiplres@cnMinLevelValF       =  signif          ; min level
  stiplres@cnMaxLevelValF       =  1.1             ; max level
  stiplres@cnLevelSpacingF      =  stiplres@cnMaxLevelValF - stiplres@cnMinLevelValF

  stiplres@cnFillColors         = (/3,3,3/)
  stiplres@cnMonoFillPattern    = False            ; want multiple patterns 
  stiplres@cnFillPatterns       = (/-1,3,-1/)     ; the patterns
  stiplres@cnMonoFillScale      = False            ; want different densities 
  stiplres@cnFillScales         = (/1.,0.5,1./)    ; the densities
  stiplres@cnFillDotSizeF       = 0.0006
  stiplres@lbLabelBarOn         = False

  ;-----------------------------
  ; Panel lay out
  ;-----------------------------
  col_lbl = (/"Revised coupling","Difference (revised - original)","Relative difference"/)
  row_lbl = species_longname

  ;-----------------------------
  ; open data file
  ;-----------------------------
  dataFileName = casename(0)+".cam.h0."+yrS+"."+dataFileTime+".180x360.nc"
  dataFile = addfile(dataFilePath(0)+dataFileName,"r")

  lat = dataFile->lat
  lon = dataFile->lon
  lev = dataFile->lev
  lev@long_name = "Norminal Pressure (hPa)"

  nlev = dimsizes(lev)
  nlat = dimsizes(lat)
  nlon = dimsizes(lon)

  ntime = dimsizes(dataFile->time)
  if (ntime.ne.1) then
     print("ntime in file /= 1")
     print("Aborting...")
     exit
  else
     itime = 0
  end if

  delete(dataFile)
  delete(dataFileName)

;=======================================
; Loop through different species
;=======================================
do irow = 0,nrow-1

  isp = irow

  print("")
  print("Processing species "+species_longname(isp))
  print("")

  ;----------------------------
  ; History variables to read 
  ;----------------------------
  if (species_category.eq.3) then
     load "./utils/mam4_species_and_modes_ax_and_cx.ncl" ; variable QoI_name becomes a list of tracer names for this species
  else
     load "./utils/mam4_species_and_modes.ncl"           ; variable QoI_name becomes a list of tracer names for this species
  end if
  nqoi = dimsizes(QoI_name)

  ;----------------------------
  ; Contour levels
  ;----------------------------

   if (.not.cvd_friendly) then
        qLevels_ref = (/1e-10,2e-10,5e-10,1e-9,2e-9,5e-9,1e-8,2e-8,5e-8,1e-7/) 
    qDiffLevels_ref = (/-3e-8,-1e-8,-3e-9,-1e-9,-3e-10,0,3e-10,1e-9,3e-9,1e-8,3e-8/)
      reldiffLevels = (/-1.5,-1,-0.5,-0.3,-0.1,0,0.1,0.3,0.5,1,1.5/)
   else
        qLevels_ref = (/3e-10,1e-9,3e-9,1e-8,3e-8,1e-7,3e-7/) 
    qDiffLevels_ref = (/-3e-8,-1e-8,-3e-9,-1e-9,-3e-10,0,3e-10,1e-9,3e-9,1e-8,3e-8/)
      reldiffLevels = (/-1.5,-1,-0.5,-0.3,-0.1,0,0.1,0.3,0.5,1,1.5/)
   end if

  if (any( species(irow).eq.(/"dst","ncl"/) )) then  
        qLevels =     qLevels_ref 
    qDiffLevels = qDiffLevels_ref 
  end if
  if (any( species(irow).eq.(/"mom","bc","pom","SOAG","SO2","H2O2"/) )) then  
        qLevels =     qLevels_ref * 1e-2
    qDiffLevels = qDiffLevels_ref * 1e-2
  end if
  if (any( species(irow).eq.(/"so4","soa"/))) then  
        qLevels =     qLevels_ref * 1e-1
    qDiffLevels = qDiffLevels_ref * 1e-1
  end if
  if (any( species(irow).eq.(/"DMS"/) )) then  
        qLevels =     qLevels_ref * 1e-3
    qDiffLevels = qDiffLevels_ref * 1e-3
  end if
  if (any( species(irow).eq.(/"H2SO4"/) )) then  
        qLevels =     qLevels_ref * 1e-5
    qDiffLevels = qDiffLevels_ref * 1e-5
  end if
  if (any( species(irow).eq.(/"O3"/) )) then  
        qLevels =     qLevels_ref * 100.
    qDiffLevels = qDiffLevels_ref * 100.
  end if

  ;-----------------------------------------------
  ; Read 3d mixing ratio fields from yearly files
  ;-----------------------------------------------
  inc_fac = 1.

  var3d_0 = new( (/nlev,nlat,nlon,nyr/),"float" )

  var3d_0!0 = "lev"
  var3d_0!1 = "lat"
  var3d_0!2 = "lon"
  var3d_0!3 = "time"
  var3d_0&lev = lev
  var3d_0&lat = lat
  var3d_0&lon = lon
    
  var3d_0 = 0.
  var3d_1 = var3d_0

  do iyr = 0,nyr-1

     dataFileName = casename+".cam.h0."+(iyr+yrS)+"."+dataFileTime+".180x360.nc"
     dataFile = addfiles(dataFilePath+dataFileName,"r")

     do iqoi = 0,nqoi-1
       varname = QoI_name(iqoi)
       var3d_0(:,:,:,iyr) = var3d_0(:,:,:,iyr) + dataFile[0]->$varname$(itime,:,:,:)*inc_fac
       var3d_1(:,:,:,iyr) = var3d_1(:,:,:,iyr) + dataFile[1]->$varname$(itime,:,:,:)*inc_fac
     end do

     print("  Done reading data from file "+dataFileName)

  end do

  res@tmYRMode = "Automatic"
  stiplres@tmYRMode = res@tmYRMode

  res@trYLog   = False
  res@trYMinF  = 300. 
  res@trYMaxF  = 1000.

 ;res@tiXAxisString = "Latitude"
  ;----------------------------------------
  ; Calculate zonal mean
  ;----------------------------------------
  var2d_0 = dim_avg_n_Wrap(var3d_0,2)  ; var3d_0's dimensions: lev,lat,lon,time
  var2d_1 = dim_avg_n_Wrap(var3d_1,2)  ; var3d_1's dimensions: lev,lat,lon,time

  ; Perform a Kolmogorov-Smirnov two-sample test the rightmost dimension (time)
  ; var2d_0's and var2d_1's dimensions: lev,lat,time

  probKS = var2d_0(:,:,0)
  probKS = (/kolsm2_n( var2d_0, var2d_1, 2 )/)

  ; Multi-year mean for plotting
  ; var2d_0's and var2d_1's dimensions: lev,lat,time

  var2d_0_mean = dim_avg_Wrap(var2d_0)
  var2d_1_mean = dim_avg_Wrap(var2d_1)

  ;===============================
  ; Plot results from case 2
  ;===============================
   res@cnLinesOn            = False
   res@gsnLeftString        = "~F22~"+row_lbl(irow)

   icol = 0
   res@cnFillPalette        = state_colormap
   res@cnLevelSelectionMode = "ExplicitLevels"   ; set explicit contour levels
   res@cnLevels             = qLevels 
   res@gsnRightString       = col_lbl(icol) 
   plot(irow,icol) = gsn_csm_pres_hgt(wks,var2d_1_mean,res)    
   delete(res@cnLevels)

  ;===============================
  ; Differences
  ;===============================
   res@cnFillPalette  = diff_colormap 
   res@cnLinesOn      = True

  ;-------------------------
  ; calculate and plot diff
  ;-------------------------
   icol = 1
   diff = var2d_1_mean
   diff = var2d_1_mean - var2d_0_mean

   res@cnLevelSelectionMode = "ExplicitLevels"   ; set explicit contour levels
   res@cnLevels             = qDiffLevels
   res@gsnRightString       = col_lbl(icol)

   plot(irow,icol) = gsn_csm_pres_hgt(wks,diff,res)
    dum(irow,icol) = gsn_csm_pres_hgt(wks,probKS,stiplres)
   overlay(plot(irow,icol),dum(irow,icol))

     ; Use dashed lines for contours of negative values, for being CVD-friendly

     cnDashPatterns = toint(where(res@cnLevels.lt.0, 2,  0  ))
     cnThicknesses  =       where(res@cnLevels.lt.0, 1., 0.1)
     set_contour_line_patterns   (plot(irow,icol),res@cnLevels,cnDashPatterns)
     set_contour_line_thicknesses(plot(irow,icol),res@cnLevels,cnThicknesses)

     delete([/cnDashPatterns,cnThicknesses/])
     delete(res@cnLevels)

  ;----------------------------------
  ; calculate and plot relative diff
  ;----------------------------------
   icol = 2 
  
   denom = where(abs(var2d_0_mean).ge.1E-36, var2d_0_mean, FillValue)
   denom@_FillValue = FillValue 
   reldiff = var2d_1_mean
   reldiff = (var2d_1_mean - var2d_0_mean)/denom

   res@cnLevelSelectionMode = "ExplicitLevels"   ; set explicit contour levels
   res@cnLevels             = reldiffLevels
   res@gsnRightString       = col_lbl(icol) 
   plot(irow,icol) = gsn_csm_pres_hgt(wks,reldiff,res)
    dum(irow,icol) = gsn_csm_pres_hgt(wks,probKS,stiplres)
   overlay(plot(irow,icol),dum(irow,icol))

     ; Use dashed lines for contours of negative values, for being CVD-friendly

     cnDashPatterns = toint(where(res@cnLevels.lt.0, 2,  0  ))
     cnThicknesses  =       where(res@cnLevels.lt.0, 1., 0.1)
     set_contour_line_patterns   (plot(irow,icol),res@cnLevels,cnDashPatterns)
     set_contour_line_thicknesses(plot(irow,icol),res@cnLevels,cnThicknesses)

     delete([/cnDashPatterns,cnThicknesses/])
     delete(res@cnLevels)

   delete([/var2d_1,var2d_0,diff,denom,reldiff/])
   delete([/var2d_1_mean,var2d_0_mean,probKS/])

delete(QoI_name)
end do

  resP = True
  resP@gsnMaximize = True
  resP@gsnPanelYWhiteSpacePercent = 5
  resP@gsnPanelXWhiteSpacePercent = 3 
  gsn_panel(wks,ndtooned(plot),(/nrow,ncol/),resP)

  print("plot file is "+plotPath+plotFile+"."+plotFmt)
