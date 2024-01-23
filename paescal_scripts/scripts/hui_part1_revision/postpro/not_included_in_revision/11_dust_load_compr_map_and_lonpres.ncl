;---------------------------------------------------------------------------------------
; This script compares the dust burden/mixing ratio from two simulations.
;  - On page 1, the accumulation mode and coarse mode are added up.
;  - If l_debug = True, two additional pages are plotted, which shows the two modes separately.
;  - Each page consists of 
;     - 3 columns: burden/mixing ration in the test simulation, diff. wrt to the control run, and rel. diff. 
;     - 4 rows: (1) map plots of vertical integrals,
;               (2) pres-lon cross section of 0-60N,
;               (3) pres-lon cross section of 40S-10S,
;               (4) pres-lat cross section of zonal mean.
;---------------------------------------------------------------------------------------
; History:
;  - Initial version by Hui Wan, PNNL, spring 2022.
;  - Revision by Hui Wan, PNNL, 2022-08, adding separate plots of the a1 and a3 modes.
;  - Revision by Hui Wan, PNNL, 2023-11, make the plots CVD-friendly.
;  - Revision by Hui Wan, PNNL, 2023-11, mask out statitistically insignificant results.
;---------------------------------------------------------------------------------------
  load "./utils/contour_line_utils.ncl"

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

  dataFilePath = "/pscratch/sd/h/huiwan/cflx/v1_papers/v1_runs/"+casename+"/climo/"

  species  = "dst"
  QoI_name = "cnd01_"+(/"dst_a1","dst_a3"/)
  nqoi     = dimsizes(QoI_name)
  chkpt    = "AERWETRM"

  if (l_debug) then
     pageTitle = "~F22~"+(/"dst_a1 + dst_a3", "dst_a1", "dst_a3"/)+"~F~"
  end if

  ;------------------
  ; contour levels
  ;------------------
    
  ;     cnLevels = (/5e-7,1e-6,2e-6,5e-6,1e-5,2e-5,5e-5,1e-4,2e-4,5e-4,1e-3/) 
  ;BdnDiffLevels = (/-5e-4,-2e-4,-1e-4,-5e-5,-2e-5,-1e-5,1e-5,2e-5, 5e-5, 1e-4,2e-4,5e-4/)

        cnLevels = (/3e-6,1e-5,3e-5,1e-4,3e-4,1e-3/)      ; fewer colors for CVD-friendly colormaps
   BdnDiffLevels = (/-3e-4,-1e-4,-3e-5,0,3e-5,1e-4,3e-4/)

  ;      qLevels = (/1e-10,2e-10,5e-10,1e-9,2e-9,5e-9,1e-8,2e-8,5e-8,1e-7/)
  ;  qDiffLevels = (/-1e-8,-5e-9,-1e-9,-5e-10,-1e-10,1e-10,5e-10,1e-9,5e-9,1e-8/)

      qLevels = (/3e-10,1e-9,3e-9,1e-8,3e-8,1e-7,3e-7/)   ; fewer colors for CVD-friendly colormaps
  qDiffLevels = (/-3e-8,-1e-8,-3e-9,-1e-9,-3e-10,0,3e-10,1e-9,3e-9,1e-8,3e-8/)

 ;relDiffLevels = (/-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8/)
  relDiffLevels = (/-0.9,-0.6,-0.3,0,0.3,0.6,0.9/)

  ;------------------
  ; color maps
  ;------------------
  cmap = read_colormap_file("batlowW")
  state_colormap = cmap(::-1,:)
  delete(cmap)
 
  cmap = read_colormap_file("bam_with_white")
  diff_colormap  = cmap(::-1,:)     ; reverse color map
  delete(cmap)

  ;-----------------------------
  ; for plotting
  ;-----------------------------
  nrow = 4  ; (1) burden map, (2)-(3) meridional avg, (4) zonal avg
  ncol = 3  ; (1) revised model, (2) revised - original, (3) relative diff.

  plot = new( (/nrow,ncol/),"graphic" )
  dum  = new( (/nrow,ncol/),"graphic" )

  plotFileTimeTag = yrRng+"_"+dataFileTime

  plotPath = "/global/cfs/projectdirs/m4359/www/huiwan/cflx/2023_rerun_for_gmd/"
  plotFile = plotPath+species+"_state_cmpr_"+plotFileTimeTag
  plotFmt  = "pdf"

  wks = gsn_open_wks(plotFmt,plotFile)
  system("chmod a+r "+plotFile+"."+plotFmt)
  gsn_define_colormap(wks,"GMT_gray")


  ;------------------------------------------
  ; For color-shaded contours
  ;------------------------------------------
  res = True
  res@gsnFrame = False
  res@gsnDraw  = False

  res@cnFillOn       = True
  res@cnLinesOn      = True

  res@tmXTOn = False
  res@tmYROn = False

  FontHeightF = 0.02 *0.36
  res@gsnStringFontHeightF = FontHeightF *1.0
  res@tiXAxisFontHeightF   = FontHeightF *0.9 
  res@tiYAxisFontHeightF   = FontHeightF *0.9

  res@lbOrientation        = "vertical"
  res@lbLabelFontHeightF   = FontHeightF *0.8

 ;res@tmXMajorGrid = True
 ;res@tmYMajorGrid = True
 ;res@tmXMajorGridThicknessF = 0.5
 ;res@tmYMajorGridThicknessF = 0.5

  ;-----------------------------------------------------
  ; For stippling that masks out insignificant results
  ;-----------------------------------------------------
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
  stiplres@cnFillScales         = (/1.,0.2,1./)    ; the densities
  stiplres@cnFillDotSizeF       = 0.0006
  stiplres@lbLabelBarOn         = False 

  ;-----------------------------
  ; Panel lay out
  ;-----------------------------

  vpWthF_map = 0.23               ; for global maps
  vpHgtF_map = vpWthF_map *0.5    ; for cross-sections

  vpWthF_crs = 0.23               ; for global maps
  vpHgtF_crs = vpWthF_crs *0.5    ; for cross-sections

  ; Set panel positions

  vpXF_map = (/0.05,0.38,0.71/)   ; for global maps
  vpXF_crs = (/0.05,0.38,0.71/)   ; for cross-sections

  vpYF_map = 0.95                 ; for global maps
  vpYF_crs = 0.95                 ; for cross-sections

  vpYF_fac = 0.65

  ; Set panel labels

  col_idx = (/"(a","(b","(c"/)
  col_lbl = (/"Revised coupling","Revised - original","Relative diff."/)
  row_lbl = (/"Column integral","0-60N avg.","40S-10S avg","Zonal avg"/)

  ; latitude bands - these should match the panel titles set above

  latmin = (/ 0.,-40./)
  latmax = (/60.,-10./)

  ;--------------------------------------
  ; Open a data file and get grid info
  ;--------------------------------------
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

  ;===============================
  if (l_debug) then
     npagebnd = dimsizes(QoI_name)
  else
     npagebnd = 0
  end if

  ;===============================
  do ipage = 0,npagebnd

     print("")
     print("Processing data for plot page "+ipage)

     ;RRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR
     ; Row 0: burden
     ;RRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR
     ; Read in the burdens calculated by CondiDiag, plot them as well as the burden differences

     irow = 0

       v_suff = "_v"
     inc_suff = ""
     inc_fac  = 1.

     ;-------------------------------------
     ; Read 2D data of individual years
     ;-------------------------------------
     if (ipage.eq.0) then  ; add all modes together and plot the sum

        iqoiS = 0
        iqoiE = dimsizes(QoI_name) -1 

     else ; read and plot a single mode

        iqoiS = ipage -1
        iqoiE = ipage -1

     end if
    
     ; Initialize arrays for summation
   
     var2d_1 = new( (/nlat,nlon,nyr/),"float" )
     var2d_1!0 = "lat"
     var2d_1!1 = "lon"
     var2d_1!2 = "time"
     var2d_1&lat = lat
     var2d_1&lon = lon
     var2d_1 = 0.

     var2d_0 = var2d_1

     ; Read and process yearly data

     do iyr = 0,nyr-1

       dataFileName = casename+".cam.h0."+(iyr+yrS)+"."+dataFileTime+".180x360.nc"
       dataFile = addfiles(dataFilePath+dataFileName,"r")

       do iqoi = iqoiS,iqoiE
          varname = QoI_name(iqoi)+v_suff+"_"+ chkpt+inc_suff
          var2d_0(:,:,iyr) = var2d_0(:,:,iyr) + dataFile[0]->$varname$(itime,:,:)*inc_fac
          var2d_1(:,:,iyr) = var2d_1(:,:,iyr) + dataFile[1]->$varname$(itime,:,:)*inc_fac
       end do 

       print("  Done reading data from file "+dataFileName)
       delete([/dataFileName,dataFile/])
     end do

     print(" Done reading data")
     ;------------------------------------------------
     ; Perform a Kolmogorov-Smirnov two-sample test 
     ;------------------------------------------------
     probKS = var2d_0(:,:,0)
     probKS = (/kolsm2_n( var2d_0, var2d_1, 2 )/)  ; KS test is applied to the rightmost dimension

     print(" Done with KS test")

     ;----------------------------------
     ; Get multi-year mean for plotting
     ;----------------------------------
     var2d_0_mean = dim_avg_Wrap(var2d_0)
     var2d_1_mean = dim_avg_Wrap(var2d_1)

     print("Calculations done for column integrals")

     ;================================================
     ; Common resources for map plots
     ;================================================
     res@gsnLeftString  = "~F22~"+row_lbl(irow)
     res@mpProjection = "CylindricalEquidistant"

     res@vpYF      = vpYF_map
     res@vpWidthF  = vpWthF_map 
     res@vpHeightF = vpHgtF_map 

     res@mpDataBaseVersion    = "LowRes"       ; better map outlines
     res@mpGridAndLimbOn      = True           ; turn on limb and grid lines
     res@mpGridLineColor      = -1             ; ...but don't draw grid lines...
     res@mpPerimOn            = False          ; turn off map perimeter
     res@mpFillDrawOrder      = "PostDraw"     ; draw map fill last
     res@mpGeophysicalLineColor = 2 ; "gray60"     ; light gray. Larger numbers are lighter.
     res@mpFillOn               = False

     res@gsnLeftString        = "~F22~"+row_lbl(irow)

     ;--------------------
     ; burden map
     ;--------------------

     res@cnLinesOn      = False
     res@cnLineLabelsOn = False

     icol = 0
     res@cnLevelSelectionMode = "ExplicitLevels"   ; set explicit contour levels
     res@cnLevels             = cnLevels 
     res@lbLabelStrings       = sprintf("%3.1e",res@cnLevels)
     res@cnFillPalette        = state_colormap 
     res@vpXF                 = vpXF_map(icol) 
     res@gsnRightString       = col_lbl(icol) 

     plot(irow,icol) = gsn_csm_contour_map(wks,var2d_1_mean,res)    

     delete(res@cnLevels)
     delete(res@lbLabelStrings)

     res@cnLinesOn      = True
     res@cnLineLabelsOn = False

     ;-------------------------
     ; burden difference map
     ;------------------------
     icol = 1

     diff = var2d_1_mean
     diff = var2d_1_mean-var2d_0_mean

     res@cnLevelSelectionMode = "ExplicitLevels"   ; set explicit contour levels
     res@cnLevels             = BdnDiffLevels 
     res@lbLabelStrings       = sprintf("%3.1e",res@cnLevels)
     res@cnFillPalette        = diff_colormap 
     res@vpXF                 = vpXF_map(icol) 
     res@gsnRightString       = col_lbl(icol) 

     plot(irow,icol) = gsn_csm_contour_map(wks,diff,res)
      dum(irow,icol) = gsn_csm_contour(wks,probKS,stiplres)
     overlay(plot(irow,icol),dum(irow,icol))

     ; use dashed lines for contours of negative values, for being CVD-friendly

     cnDashPatterns = toint(where(res@cnLevels.lt.0, 2,  0  ))
     cnThicknesses  =       where(res@cnLevels.lt.0, 1., 0.1)
     set_contour_line_patterns   (plot(irow,icol),res@cnLevels,cnDashPatterns)
     set_contour_line_thicknesses(plot(irow,icol),res@cnLevels,cnThicknesses)

     delete(res@cnLevels)
     delete(res@lbLabelStrings)
     delete([/cnDashPatterns,cnThicknesses/])

     ;--------------------------------
     ; burden relative difference map
     ;--------------------------------
     icol = 2
     denom = where(abs(var2d_0_mean).ge.1E-36, var2d_0_mean, -999.)
     denom@_FillValue = -999.

     reldiff = var2d_1_mean 
     reldiff = ( var2d_1_mean - var2d_0_mean )/denom

     res@cnLevelSelectionMode = "ExplicitLevels"   ; set explicit contour levels
     res@cnLevels             = relDiffLevels
     res@vpXF                 = vpXF_map(icol) 
     res@gsnRightString       = col_lbl(icol) 

     plot(irow,icol) = gsn_csm_contour_map(wks,reldiff,res)
      dum(irow,icol) = gsn_csm_contour(wks,probKS,stiplres)
     overlay(plot(irow,icol),dum(irow,icol))

     ; Use dashed lines for contours of negative values, for being CVD-friendly

     cnDashPatterns = toint(where(res@cnLevels.lt.0, 2,  0  ))
     cnThicknesses  =       where(res@cnLevels.lt.0, 1., 0.1)
     set_contour_line_patterns   (plot(irow,icol),res@cnLevels,cnDashPatterns)
     set_contour_line_thicknesses(plot(irow,icol),res@cnLevels,cnThicknesses)

     delete([/cnDashPatterns,cnThicknesses/])
     delete(res@cnLevels)

     ;--------------------------------
     ; Done with map plots. Clean up
     ;--------------------------------
     delete([/var2d_1,var2d_0,diff,denom,reldiff/])
     delete([/var2d_1_mean,var2d_0_mean,probKS/])
     delete(res@mpProjection         )
     delete(res@mpDataBaseVersion    )
     delete(res@mpGridAndLimbOn      )
     delete(res@mpGridLineColor      )
     delete(res@mpPerimOn            )
     delete(res@mpFillDrawOrder      )
     delete(res@mpGeophysicalLineColor)
     delete(res@mpFillOn)

     print("Plotting done for column integrals")
     print("")

     ;RRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR
     ; Rows 1-3: cross sections.
     ;RRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR
     ; Initialize arrays to be for summation
   
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

     ; read and process data

     v_suff = ""
     inc_fac = 1.

     do iyr = 0,nyr-1

       dataFileName = casename+".cam.h0."+(iyr+yrS)+"."+dataFileTime+".180x360.nc"
       dataFile = addfiles(dataFilePath+dataFileName,"r")

       do iqoi = iqoiS,iqoiE
          varname = QoI_name(iqoi)+v_suff+"_"+ chkpt+inc_suff
          var3d_0(:,:,:,iyr) = var3d_0(:,:,:,iyr) + dataFile[0]->$varname$(itime,:,:,:)*inc_fac
          var3d_1(:,:,:,iyr) = var3d_1(:,:,:,iyr) + dataFile[1]->$varname$(itime,:,:,:)*inc_fac
       end do 

       print("  Done reading data from file "+dataFileName)
       delete([/dataFileName,dataFile/])
     end do

     ; Swap the eastern and western hemispheres so that lon-pressure cross-sections
     ; will have lon = 0 in the middle

     lontmp = lon
     lontmp(0:nlon/2-1) = lon(nlon/2:) 
     lontmp(nlon/2:)    = lon(0:nlon/2-1)
     lontmp = where(lon.gt.180., lon-360., lon)

     var3d_0&lon = lontmp
     var3d_1&lon = lontmp

     tmp = var3d_0
     var3d_0(:,:,0:nlon/2-1,:) = tmp(:,:,nlon/2:,   :)
     var3d_0(:,:,nlon/2:,   :) = tmp(:,:,0:nlon/2-1,:)

     tmp = var3d_1
     var3d_1(:,:,0:nlon/2-1,:) = tmp(:,:,nlon/2:,   :)
     var3d_1(:,:,nlon/2:,   :) = tmp(:,:,0:nlon/2-1,:)
     delete(tmp)

     ;RRRRRRRRRRRRRRRRRRRRRRRRRR
     ; Zonal mean
     ;RRRRRRRRRRRRRRRRRRRRRRRRRR
     irow = 3

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

     ; Plot settings

     stiplres@tmYRMode = "Automatic"
     res@tmYRMode = "Automatic"
     res@trYLog   = False
     res@trYMinF  = 100. 
     res@trYMaxF  = 1000.

     res@trXMinF = -90.
     res@trXMaxF =  90.

     res@vpWidthF  = vpWthF_crs 
     res@vpHeightF = vpHgtF_crs 

     res@tiXAxisString = "Latitude"
     res@vpYF = vpYF_crs - irow * vpHgtF_crs/vpYF_fac

     res@gsnLeftString  = "~F22~"+row_lbl(irow)

     ;---------------------------------------
     ; Plot cross section of the field itself
     ;---------------------------------------
     res@cnLinesOn      = False
     res@cnLineLabelsOn = False

     icol = 0
     res@cnFillPalette        = state_colormap
     res@cnLevelSelectionMode = "ExplicitLevels"   ; set explicit contour levels
     res@cnLevels             = qLevels 
     res@vpXF                 = vpXF_crs(icol) 
     res@gsnRightString       = col_lbl(icol)

     plot(irow,icol) = gsn_csm_pres_hgt(wks,var2d_1_mean,res)
     delete(res@cnLevels)

     ;----------------------------------
     ; Calculate and plot diff
     ;----------------------------------

     res@cnLinesOn      = True
     res@cnLineLabelsOn = False
     res@cnFillPalette  = diff_colormap 

     icol = 1
     diff = var2d_1_mean
     diff = var2d_1_mean-var2d_0_mean
    ;diff = where( probKS.lt.signif, var2d_1_mean-var2d_0_mean, FillValue )
     diff@_FillValue = FillValue

     res@cnLevelSelectionMode = "ExplicitLevels"   ; set explicit contour levels
     res@cnLevels             = qDiffLevels 
     res@vpXF                 = vpXF_crs(icol) 
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
     ; Calculate and plot relative diff
     ;----------------------------------
     icol = 2 
  
     denom = where(abs(var2d_0_mean).ge.1E-36, var2d_0_mean, -999.)
     denom@_FillValue = -999.

     reldiff = var2d_1_mean 
     reldiff = (var2d_1_mean - var2d_0_mean)/denom
    ;reldiff = where( probKS.lt.signif, (var2d_1_mean - var2d_0_mean)/denom, FillValue )
     reldiff@_FillValue = FillValue

     res@cnLevelSelectionMode = "ExplicitLevels"   ; set explicit contour levels
     res@cnLevels             = relDiffLevels
     res@vpXF                 = vpXF_crs(icol) 
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

     ;RRRRRRRRRRRRRRRRRRRRRRRRRR
     ; Meridional mean
     ;RRRRRRRRRRRRRRRRRRRRRRRRRR
     res@tiXAxisString = "Longitude"
     res@trXMinF = -180.
     res@trXMaxF =  180.

     do irow = 1,2

        var2d_0 = dim_avg_n_Wrap(var3d_0(:,{latmin(irow-1):latmax(irow-1)},:,:),1)  ; var3d_0's dimensions: lev,lat,lon,time 
        var2d_1 = dim_avg_n_Wrap(var3d_1(:,{latmin(irow-1):latmax(irow-1)},:,:),1)  ; var3d_1's dimensions: lev,lat,lon,time 

        ; Perform a Kolmogorov-Smirnov two-sample test the rightmost dimension (time)
        ; var2d_0's and var2d_1's dimensions: lev,lon,time

        probKS = var2d_0(:,:,0)
        probKS = (/kolsm2_n( var2d_0, var2d_1, 2 )/)

        ; Multi-year mean for plotting
        ; var2d_0's and var2d_1's dimensions: lev,lon,time

        var2d_0_mean = dim_avg_Wrap( var2d_0)
        var2d_1_mean = dim_avg_Wrap( var2d_1)


        res@vpYF = vpYF_crs - irow * vpHgtF_crs/vpYF_fac
        res@gsnLeftString = "~F22~"+row_lbl(irow)

        ;---------------------------------------
        ; Plot cross-section of the field itself
        ;---------------------------------------
        res@cnLinesOn      = False
        res@cnLineLabelsOn = False

        icol = 0
        res@cnFillPalette        = state_colormap
        res@cnLevelSelectionMode = "ExplicitLevels"   ; set explicit contour levels
        res@cnLevels             = qLevels 
        res@vpXF                 = vpXF_crs(icol)
        res@gsnRightString       = col_lbl(icol)

        plot(irow,icol) = gsn_csm_pres_hgt(wks,var2d_1_mean,res)    
        delete(res@cnLevels)

        ;-------------------------
        ; calculate and plot diff
        ;-------------------------
        res@cnLinesOn      = True
        res@cnLineLabelsOn = False


        icol = 1
        res@cnFillPalette  = diff_colormap

        diff = var2d_1_mean
        diff = var2d_1_mean-var2d_0_mean
       ;diff = where( probKS.lt.signif, var2d_1_mean-var2d_0_mean, FillValue )
        diff@_FillValue = FillValue

        res@cnLevelSelectionMode = "ExplicitLevels"   ; set explicit contour levels
        res@cnLevels             = qDiffLevels
        res@gsnRightString       = col_lbl(icol)
        res@vpXF                 = vpXF_crs(icol)

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

        ;-------------------------------------
        ; Calculate and plot relative diff
        ;-------------------------------------

        icol = 2
        res@cnFillPalette  = diff_colormap

        res@cnLevelSelectionMode = "ExplicitLevels"   ; set explicit contour levels
        res@cnLevels  = relDiffLevels 

        denom = where(abs(var2d_0_mean).ge.1E-36, var2d_0_mean, -999.)
        denom@_FillValue = -999.
        reldiff = var2d_1_mean 
        reldiff = (var2d_1_mean-var2d_0_mean)/denom
       ;reldiff = where( probKS.lt.signif, (var2d_1_mean-var2d_0_mean)/denom, FillValue )
        reldiff@_FillValue = FillValue

        res@gsnRightString       = col_lbl(icol)
        res@vpXF                 = vpXF_crs(icol)
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

     end do
     delete(res@tiXAxisString)
     delete(res@trXMinF)
     delete(res@trXMaxF)
     delete(res@trYMinF)
     delete(res@trYMaxF)

     if (l_debug) then ; add a page title
        txres = True
        txres@txFontHeightF = FontHeightF*1.2
        gsn_text_ndc(wks,pageTitle(ipage),0.5,vpYF_map+0.03,txres)
     end if
 
     ;----------------------
     ; Draw panels on page
     ;----------------------
     psres = True
     maximize_output(wks,psres)  ; calls draw and frame for you 

  end do ; page loop
  ;===========================================================================
