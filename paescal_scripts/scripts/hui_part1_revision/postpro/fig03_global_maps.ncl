  ;**********************************************************************************************************************
  ; Plot global maps of emissions and burden
  ;**********************************************************************************************************************
  ; History files with global output

  dataFileName = casename+".cam.h0."+yrRng+"."+dataFileTime+".nc"
  dataFile = addfiles( dataFilePath+dataFileName,"r")

  plotFileBasename = species+"_global_maps"

  if (l_debug) then
     caseloopbnd = ncase-1    ; plot all cases for examination
  else
     caseloopbnd = 0          ; plot only the first case for paper 
  end if

  do icase = 0,caseloopbnd ;#######################################

     plotFileSuffix = "_"+plotFileTimeTag+"_"+casename(icase)
     print("Plotting global maps for "+plotFileSuffix)

     ; Create a plot file

     plotFile = plotPath+plotFileBasename+plotFileSuffix
     plotFmt  = "pdf"
     wks      = gsn_open_wks(plotFmt,plotFile)
     system("chmod a+r "+plotFile+"."+plotFmt)
     gsn_define_colormap(wks,"BatlowW")
     gsn_reverse_colormap(wks)

     ; We will only plot 2 rows and 1 column, but in order to have similar font sizes in 
     ; the map plot and profile plot, we declare the variable plot to have the same shape 
     ; as the variable xyplot. 

     global_map = new( 2,"graphic" )
     dum4box    = new( (/2,nregion,4/),"graphic")

     ; Read in lat/lon; create tmp array (as we will need to sum over different QoIs - modes in this case)

     cnres@sfXArray = dataFile[icase]->lon
     cnres@sfYArray = dataFile[icase]->lat

     ncol  = dimsizes(cnres@sfXArray)
     fld2d = new(ncol,"float")

     ;---------------------------------------------------------------
     ; First map plot: emissions with boxes indicating focus regions
     ;---------------------------------------------------------------
     cnres@tiMainString = "(a) Dust emissions"
     cnres@gsnCenterString = "Mass flux, unit: kg m~S~-2~N~ s~S~-1~N~"+"~C~"
     cnres@cnLevelSelectionMode = "ExplicitLevels"   ; set explicit contour levels
    ;cnres@cnLevels = (/5e-11, 1e-10, 5e-10, 1e-9, 5e-9, 1e-8, 5e-8 /)   ; more colors
     cnres@cnLevels = (/1e-12,1e-11,1e-10,1e-9,1e-8,1e-7/)          ; fewer colors when using CVD-friendly colormaps

     ; Read data and make a map plot

     ichkpt = ichkpt_emis  ; emissions
     inc    = "_inc"       ; inc/dtime gives tendency
     scalef = 1./dtime

     fld2d  = 0.
     do iqoi = 0,nqoi-1
        varname = QoI_name(iqoi)+"_v_"+ chkpt(icase,ichkpt)+inc
        fld2d = fld2d + dataFile[icase]->$varname$(itime,:) * scalef 
     end do
     global_map(0) = gsn_csm_contour_map(wks,fld2d,cnres)

     ; Add boxes to show the regions

     respoly                  = True                      ; polyline mods desired
     respoly@gsLineThicknessF = 2.0
     box_LineDashPattern      = (/0,0,11/)
     do ir = 0,nregion-1
        xpts = (/lonmin(ir),lonmax(ir),lonmax(ir),lonmin(ir),lonmin(ir)/)
        ypts = (/latmin(ir),latmin(ir),latmax(ir),latmax(ir),latmin(ir)/)
        respoly@gsLineColor       = 1 ; black
        respoly@gsLineDashPattern = box_LineDashPattern(ir)       
        do ii = 0,3
           dum4box(0,ir,ii) = gsn_add_polyline(wks,global_map(0),xpts(ii:ii+1),ypts(ii:ii+1),respoly)
        end do
     end do

     delete(cnres@cnLevels)
     ;-----------------------------------------
     ; Second map plot: burden (and no boxex)
     ;-----------------------------------------
     cnres@tiMainString = "(b) Dust burden"
     cnres@gsnCenterString = "Mass load, unit: kg m~S~-2~N~"+"~C~"
    ;cnres@cnLevels = (/5e-7,1e-6,2e-6,5e-6,1e-5,2e-5,5e-5,1e-4,2e-4,5e-4,1e-3/)  ; more colors
     cnres@cnLevels = (/1e-6,3e-6,1e-5,3e-7,1e-4,3e-4,1e-3/) ; fewer colors when using CVD-friendly colormaps

     ; Read data and make a map plot

     ichkpt = ichkpt_wtrm  ; wet removal
     inc    = ""           ; "" means state 
     scalef = 1.

     fld2d  = 0.
     do iqoi = 0,nqoi-1
        varname = QoI_name(iqoi)+"_v_"+ chkpt(icase,ichkpt)+inc
        fld2d = fld2d + dataFile[icase]->$varname$(itime,:) * scalef
     end do
     global_map(1) = gsn_csm_contour_map(wks,fld2d,cnres)

     ;-----------------------------------------
     ; Panel plot
     ;-----------------------------------------
     resP = True
     resP@gsnPanelYWhiteSpacePercent = 10 
     resP@gsnPanelXWhiteSpacePercent = 10
     gsn_panel(wks,global_map,(/1,2/),resP)

     delete(cnres@sfXArray)
     delete(cnres@sfYArray)
     delete(cnres@cnLevels)
     delete(fld2d)
     delete(resP)

  end do ;#######
  delete(dataFile)
  delete(dataFileName)

  print("")
  print("Done plotting global maps")
  print("***************************")
  print("")
