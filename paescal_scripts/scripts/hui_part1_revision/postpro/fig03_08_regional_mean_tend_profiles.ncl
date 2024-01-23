
  ; Open one data file and get grid info

  dataFileName = casename(0)+".cam.h0."+yrS+"."+dataFileTime+".nc"
  dataFile = addfile(dataFilePath(0)+dataFileName,"r")

  ;------------------------------------------
  ; Need to get nlev already now
  ;------------------------------------------
  preslev = dataFile->lev
  preslev@long_name = "Norminal pressure (hPa)"
  nlev = dimsizes(preslev)

  FillValue = -999.
  region_avg_single_yr = new( (/ncase,nregion,nchkpt,nlev,nyr/), "float" )
  region_avg_single_yr@_FillValue = FillValue

  ;---------------------------------------------------
  ; Get lon/lat/area since they don't change
  ;---------------------------------------------------
  lon  = dataFile->lon
  lat  = dataFile->lat
  area = tofloat(dataFile->area)
  ncol = dimsizes(lon)

  if (l_debug) then
      print("")
      print("lon range in file: "+min(lon)+", "+max(lon))
      print("lat range in file: "+min(lat)+", "+max(lat))
      print("")
  end if

  ;#########################################################################
  ; Read in data and calculate regional averages for each case (simulation)
  ;#########################################################################
  do icase = 0,ncase-1

     print("===============================================================")
     print(" casename = "+casename(icase))
     print("===============================================================")

     ;---------------------------------------------------------------------------------------
     ; For debugging, create a plot file for plotting the regions on a map for sanity check 
     ;---------------------------------------------------------------------------------------
     if (l_debug) then

        map_plot = new( nregion, "graphic" )

        plotFile = plotPath+"/regions_"+plotFileTimeTag+"_"+casename(icase)
        plotFmt  = "pdf"
        wks0     = gsn_open_wks(plotFmt,plotFile)
        system("chmod a+r "+plotFile+"."+plotFmt)
        gsn_define_colormap(wks0,"testcmap")

     end if ; l_debug

     ;-----------------------------------------------------------------------------
     ; Do the calculation for each year separately as we need stddev for the plots 
     ;----------------------------------------------------------------------------
     do iyr = 0,nyr-1

        ;----------------
        ; History files
        ;----------------
        yyyy = iyr + yrS
        dataFileName = casename(icase)+".cam.h0."+yyyy+"."+dataFileTime+".nc"
        dataFile = addfile( dataFilePath(icase)+dataFileName,"r")
        print("processing data from file "+dataFileName)

        ntime = dimsizes(dataFile->time)

        if (ntime.ne.1) then
           print("ntime in file /= 1")
           print("Aborting...")
           exit
        else
           itime = 0
        end if

        ;---------------------------------------------------------------------------------------------
        ; Use total emission of all modes for masking. So, read in the global fields of emisson now
        ; and calculate the sum for this species.
        ;---------------------------------------------------------------------------------------------
        ichkpt = ichkpt_emis  ; emissions
        inc    = "_inc"       ; inc/dtime gives tendency
        scalef = 1./dtime     ; convert from burden increment to emission flux

        emissions = new(ncol,"float")  ; Create tmp array, as we will need to sum over different QoIs (i.e., modes in this case)
        emissions = 0.
        do iqoi = 0,nqoi-1
           varname = QoI_name(iqoi)+"_v_"+ chkpt(icase,ichkpt)+inc
           emissions = emissions + dataFile->$varname$(itime,:) * scalef
        end do
        print("done reading emissions")

        ;-------------------------------------------------------------------------------------------------------
        ; Select grid cells in each region, then read in increments/tendencies and calculate spatial avg.
        ;-------------------------------------------------------------------------------------------------------
        emis_cutoff = 0.

        do ir = 0,nregion-1

           ; Find grid cells in the lat lon box; get their indices

           latlon_flag  = cell_is_in_latlon_box( lat,lon,latmin(ir),latmax(ir),lonmin(ir),lonmax(ir) ) ; returns a 0/1 flag

           idx_cell_in_latlon_box = ind(latlon_flag.eq.1)
           found_cells = .not.(all(ismissing(idx_cell_in_latlon_box)))

           if (l_debug) then
              print("Region "+ir+", found cells in latlon range = "+found_cells)
           end if

           ; Check whether there is emission in each grid cell

           emission_flag = where(emissions.gt.emis_cutoff, 1, 0)

           ; Select grid cells that are both of the emission type of interest and fall into the lat lon box

           cell_selected = ind( (emission_flag.eq.emis_type(ir)) .and. (latlon_flag.eq.1) )

           if (found_cells) then
             ;------------------------------------------------------------
             ; Plot values in latlon boxes for sanity check
             ;------------------------------------------------------------
             if (l_debug.and.(iyr.eq.0)) then

              cnres@sfXArray = lon(idx_cell_in_latlon_box)
              cnres@sfYArray = lat(idx_cell_in_latlon_box)
              cnres@tiMainString = "Region "+ir
              cnres@gsnCenterString = "lon + lat"+"~C~"

              var2d = lon(idx_cell_in_latlon_box)+lat(idx_cell_in_latlon_box)
              map_plot(ir) = gsn_csm_contour_map(wks0,var2d,cnres)

              delete(var2d)
              delete(cnres@sfXArray)
              delete(cnres@sfYArray)

             end if ; l_debug
             ;----------------

             ;---------------------------------------------------------------------
             ; Calculate regional averages of total tendencies of all modes (QoIs)
             ;---------------------------------------------------------------------
             ; Note that different checkpoints and levels are processed separately but using the same masking.
             ; So, essentially, we average over all grid columns region to get vertical profiles of all
             ; processes of interest.

             inc_suff = "_inc"      ; increment
             inc_fac  = 1./dtime    ; convert to tendency

             do ichkpt = 0,nchkpt-1

                ; Read in tendencies and sum up all modes

                var3d = new( (/nlev,ncol/),"float" )
                var3d = 0.

                do iqoi = 0,nqoi-1
                   varname = QoI_name(iqoi)+"_"+ chkpt(icase,ichkpt)+inc_suff
                   var3d = var3d + dataFile->$varname$(itime,:,:)*inc_fac
                end do

                ; Average over region.

                do ilev = 0,nlev-1
                   region_avg_single_yr(icase,ir,ichkpt,ilev,iyr) = sum( area(cell_selected)*var3d(ilev,cell_selected) ) \
                                                                   /sum( area(cell_selected) )
                end do

                if (l_debug) then
                   print(" ** checkpoint "+chkpt(icase,ichkpt)+inc_suff+" region avg done")
                end if
             end do  ; chkpt loop

           end if ; if found_cells

           delete([/latlon_flag,idx_cell_in_latlon_box,emission_flag,cell_selected/])
        end do ; region

     end do ; year

     ;------------------------------------------------------------
     ; Panel the map plots of latlon boxes for sanity check
     ;------------------------------------------------------------
     if (l_debug) then
        gsn_panel(wks0,ndtooned(map_plot),(/nregion,1/),False)
        delete(emissions)
     end if
     print("")

  end do ; case loop

  region_avg_mean = dim_avg   ( region_avg_single_yr )  ; multi-year mean
  region_avg_std  = dim_stddev( region_avg_single_yr )  ; std. dev. of individual years

  print("Calculation of regional averages done.")
  print("")
  print("**************************************")
  print("")


  ;###############################################################################################
  ; Plot the regional averages calculated above. Global fields are no longer needed.
  ;###############################################################################################
  ;
  ; We create one plot file for each case (simulation) that contains all profiles and regions.
  ; After that, we create a 3rd file that has different simulations in one figure for comparison.
  ;
  ; In each of these files, the profiles are plotted twice, i.e., in two rows focusing on different vertical layers.
  ;------------------------------------------------------------------------------------------------------------------

  lev_opt_label = (/"Bottom layer excluded","Lowest 10 layers"/)
  nlev_opt = dimsizes(lev_opt_label)
  xyplot = new( (/nlev_opt,nregion/),"graphic" )                 ; multi-year mean
  plydum = new( (/nlev_opt,nregion,nchkpt*2,nlev*2/),"graphic" ) ; variability

  ; Panel titles

  fig_col_title_prefix = new( (/nlev_opt,nregion/),"string" )
  fig_col_title_prefix(:,0) = region_label(0)+" mean tendencies"
  fig_col_title_prefix(:,1) = region_label(1)+" mean tendencies"
  fig_col_title_prefix(:,2) = region_label(2)+" mean tendencies"

  ; Manually set scale factors to accommodate the very different magnitudes of tendencies
  ; in different types of regions and at different altitudes.

  scale_fac = new( (/nlev_opt,nregion/),"float" )
  scale_fac(0,0) = 100.
  scale_fac(1,0) = 1.
  scale_fac(0,1) = 1000.
  scale_fac(1,1) = 100.
  scale_fac(0,2) = 100000.
  scale_fac(1,2) = 10000.

  ;------------------------------------------------------------------------------------------------------------------
  plotFileBasename = species+"_tend_profiles"

  ; Save a copy of the tendency profiles, as we will apply fill values in some panels to mask out
  ; emission-induced tendencies that are zero (to address a comment from reviewer 2 of the manuscript).

  region_avg_mean_save = region_avg_mean

  ; Tendencies are small. So scale them up for plotting so that we don't need to write x10^* for the tick marks.

  unit     = "10~S~-10~N~ kg kg~S~-1~N~ s~S~-1~N~"
  unit_fac = 1e10

  ;--------------------------------------------

  res = True
  res@gsnFrame = False
  res@gsnDraw  = False
  res@vpHeightF = 0.6

  res@tiMainFont = 22 ;helvetica-bold
  res@tiMainFontHeightF    = FontHeightF*1.1
  res@gsnStringFontHeightF = FontHeightF
  res@tiXAxisFontHeightF   = FontHeightF
  res@tiYAxisFontHeightF   = FontHeightF
  res@tmXBLabelFontHeightF = FontHeightF
  res@tmYLLabelFontHeightF = FontHeightF

  res@tmXTOn = False
  res@tmYROn = False
  res@tmYLFormat = "f"
  res@trYReverse = True

  res@gsnXRefLine = 0.
  res@gsnXRefLineDashPattern = 0
  res@gsnXRefLineColor  = "black"

  res@xyMarkLineMode = "MarkLines"

  res@pmLegendDisplayMode    = "Always"
  res@pmLegendSide           = "Top"
  res@pmLegendWidthF         = 0.08

  ;-----------------------------------------------------------------------------------------
  ;-----------------------------------------------------------------------------------------
  caseloopbnd = ncase
  do icase = 0,caseloopbnd

    l_skip = (.not.l_debug).and.(icase.eq.1)  ; skip plotting the second case unless in debug mode
    if (.not.l_skip) then

     l_compare_cases = icase.ge.ncase

     ;-------------------------------------------
     ; Create a plot file
     ;-------------------------------------------
     if (l_compare_cases) then
        print("Plotting a comparison of profiles from cases ")
        print("  "+casename)
        plotFile = plotPath+plotFileBasename+"_cmpr_"+plotFileTimeTag
     else
        print("Plotting profiles for case "+icase)
        plotFile = plotPath+plotFileBasename+"_"+plotFileTimeTag+"_"+casename(icase)
     end if

     plotFmt = "pdf"
     wks = gsn_open_wks(plotFmt,plotFile)
     system("chmod a+r "+plotFile+"."+plotFmt)
     gsn_define_colormap(wks,"StepSeq25")

     colors = (/4,14,19,25,9/)

     ; Marker style settings. Use different marker shapes and fill-vs-open to be CVD friendly.
     ; This has to be done here, after "wks = gsn_open_wks(plotFmt,plotFile)" because the new
     ; markers are associated with wks.

     wksLocal = wks
     load "./utils/mymarkers.ncl"
     delete(wksLocal)
     markers  = (/MarkerFildLeaf,   MarkerFildCircle, MarkerFildSquare, MarkerFildTriagl, MarkerFildStar/)
     markers2 = (/MarkerOpenRhombd, MarkerOpenCircle, MarkerOpenSquare, MarkerOpenTriagl, MarkerOpenStar/)

     ;----------------------------------------------------------------------
     ; Different different columns of panels show different regions.
     ; Different different rows show different altitude ranges.
     ;----------------------------------------------------------------------
     do iregion = 0,nregion-1  
     do ilev_opt = 0,nlev_opt-1

           ;---------------------------------------------------------------
           ; In the upper and lower rows, we label the y-axis differently
           ;---------------------------------------------------------------
           if (ilev_opt.eq.0) then  ; upper row: show pressure levels; do not show lowest level

              res@xyMarkerSizeF = res@vpHeightF * 0.022
              res@trYMaxF = 1020.
              res@trYMinF = 500.
              lev = preslev
              plotlevmin = 0 
              plotlevmax = nlev-2

           else ; lower row: show layer indices; include the lowest layer
 
              res@xyMarkerSizeF = res@vpHeightF * 0.035
              res@trYMaxF = nlev + 1 
              res@trYMinF = nlev - 10 
              lev = ispan(1,nlev,1)
              lev@long_name = "Model layer index"
              plotlevmin = nlev-10 
              plotlevmax = nlev-1
           end if

           ;---------------------------------------------------------------
           ; User a larger x-axis range when comparing two simulations
           ;---------------------------------------------------------------
           if (l_compare_cases) then
             res@trXMaxF =  8. 
             res@trXMinF = -8. 
           else
             res@trXMaxF =  5. 
             res@trXMinF = -6.
           end if 

           res@tmXBMode = "Manual"
           res@tmXBTickSpacingF = res@trXMaxF/2.

           ;-----------------------------------------------------------------------------------------------------
           ; This is quite tedious. Better unify the handling of the two situations (single exp vs comparison).
           ;-----------------------------------------------------------------------------------------------------
           if (l_compare_cases) then

              ; for mean profiles

              nk = dimsizes(lev(plotlevmin:plotlevmax))
              xx = new( (/2*(nchkpt-1),nk/), "float")
              xx(0:nchkpt-2,:) = region_avg_mean(0,iregion,1:,plotlevmin:plotlevmax) * scale_fac(ilev_opt,iregion) * unit_fac
              xx(nchkpt-1:, :) = region_avg_mean(1,iregion,1:,plotlevmin:plotlevmax) * scale_fac(ilev_opt,iregion) * unit_fac
              yy = lev(plotlevmin:plotlevmax)

              ; for variability ranges

              npr = 2*(nchkpt-1)   ; # of profiles to plot range for
              nyy = dimsizes(yy)

              prfPoly = new( (/npr,nyy,2/),typeof(region_avg_mean) )

              prfPoly(0:nchkpt-2,:,0) =( region_avg_mean(0,iregion,1:,plotlevmin:plotlevmax) \
                                        -region_avg_std (0,iregion,1:,plotlevmin:plotlevmax) ) \
                                       *scale_fac(ilev_opt,iregion) * unit_fac

              prfPoly(0:nchkpt-2,:,1) =( region_avg_mean(0,iregion,1:,plotlevmin:plotlevmax) \
                                        +region_avg_std (0,iregion,1:,plotlevmin:plotlevmax) ) \
                                       *scale_fac(ilev_opt,iregion) * unit_fac

              prfPoly(nchkpt-1:,:,0)  =( region_avg_mean(1,iregion,1:,plotlevmin:plotlevmax) \
                                        -region_avg_std (1,iregion,1:,plotlevmin:plotlevmax) ) \
                                       *scale_fac(ilev_opt,iregion) * unit_fac

              prfPoly(nchkpt-1:,:,1)  =( region_avg_mean(1,iregion,1:,plotlevmin:plotlevmax) \
                                        +region_avg_std (1,iregion,1:,plotlevmin:plotlevmax) ) \
                                       *scale_fac(ilev_opt,iregion) * unit_fac

              ;---

              res@xyExplicitLegendLabels = new( 2*(nchkpt-1), "string" )
              res@xyExplicitLegendLabels(:nchkpt-2)  = "  "+label(1:)+" (ori) "
              res@xyExplicitLegendLabels(nchkpt-1:) = "  "+label(1:)+" (rev) "

              res@xyMarkerColors = new( 2*(nchkpt-1),"integer" )
              res@xyMarkerColors(0:nchkpt-2) = colors(1:nchkpt-1)   ; use the same color for the same physical process
              res@xyMarkerColors(nchkpt-1:)  = colors(1:nchkpt-1)

              res@xyDashPatterns = new( 2*(nchkpt-1),"integer" )
              res@xyDashPatterns(0:nchkpt-2) = 0   ; use different dash patterns to distinguishi cases
              res@xyDashPatterns(nchkpt-1:)  = 1

              res@xyMarkers = new( 2*(nchkpt-1), "integer" )
              res@xyMarkers(0:nchkpt-2) = markers (1:nchkpt-1)   ; use different markers (filled or unfilled) to further help distinguish cases
              res@xyMarkers(nchkpt-1:)  = markers2(1:nchkpt-1) 

              res@xyMarkerThicknessF = 0.7
              res@vpHeightF = 0.6
              res@pmLegendHeightF = res@vpHeightF *0.4
              res@pmLegendOrthogonalPosF = -0.6
              res@pmLegendParallelPosF   = .22
              res@lgLabelFontHeightF     = res@tmYLLabelFontHeightF*0.6

              fig_col_label = fig_col_title_prefix
              fig_col_label = ""

           else
              region_avg_mean(icase,iregion,:,:) = region_avg_mean_save(icase,iregion,:,:)

             ; Per reviewer's request, do not plot emission-caused tendencies when they are zero
             if ( emis_type(iregion).ne.1 .or. \
                 (emis_type(iregion).eq.1).and.(ilev_opt.eq.0)) then
                region_avg_mean(icase,iregion,ichkpt_emis,:) = FillValue
             end if

              ; for mean profiles

              xx = region_avg_mean(icase,iregion,:,plotlevmin:plotlevmax) * scale_fac(ilev_opt,iregion) * unit_fac
              yy = lev(plotlevmin:plotlevmax)

              ; for variability ranges

              npr = nchkpt       ; # of profiles to plot range for
              nyy = dimsizes(yy)

              prfPoly = new( (/npr,nyy,2/),typeof(region_avg_mean) )

              prfPoly(:,:,0) =  ( region_avg_mean(icase,iregion,:,plotlevmin:plotlevmax) \
                                 -region_avg_std (icase,iregion,:,plotlevmin:plotlevmax) ) \
                               *scale_fac(ilev_opt,iregion) * unit_fac

              prfPoly(:,:,1) =  ( region_avg_mean(icase,iregion,:,plotlevmin:plotlevmax) \
                                 +region_avg_std (icase,iregion,:,plotlevmin:plotlevmax) ) \
                               *scale_fac(ilev_opt,iregion) * unit_fac

              ;--------

              res@xyExplicitLegendLabels = "  "+label+" "
              res@xyMarkerColors = colors
              res@xyDashPattern  = 0
              res@xyMarkers      = markers
              res@xyMarkerThicknessF = 1.
              res@vpHeightF      = 0.6
              res@pmLegendHeightF = res@vpHeightF *0.25 ; height of legend.
              res@pmLegendOrthogonalPosF = res@vpHeightF*(-0.7)               ; move units down
              res@pmLegendParallelPosF   = .24                 ; move units right
              res@lgLabelFontHeightF     = res@tmYLLabelFontHeightF*0.8

              fig_col_label = fig_col_title_prefix

              fig_col_label(0,0) = "(c)"
              fig_col_label(0,1) = "(d)"
              fig_col_label(0,2) = "(e)"
              fig_col_label(1,0) = "(f)"
              fig_col_label(1,1) = "(g)"
              fig_col_label(1,2) = "(h)"
           end if

           res@xyLineColors   = res@xyMarkerColors

           if (scale_fac(ilev_opt,iregion).le.1e3) then
              scale_fac_string = " x"+scale_fac(ilev_opt,iregion)+"    " 
           else
              scale_fac_string = " x"+sprintf("%.0g", scale_fac(ilev_opt,iregion) )+"     "
           end if
           res@tiMainString = fig_col_label(ilev_opt,iregion)+" "+\
                              fig_col_title_prefix(ilev_opt,iregion)+\
                              scale_fac_string

           res@gsnCenterString = lev_opt_label(ilev_opt)
           res@tiXAxisString = "Mass mixing ratio tendency ("+unit+")"

           xyplot(ilev_opt,iregion) = gsn_csm_xy(wks,xx,yy,res)

           ;-------------------------------------------
           ; Add horizontal bars to show variability
           ;-------------------------------------------
           gsres = True
           gsres@tfPolyDrawOrder  = "Postdraw" ; put this on top of markers/lines
           gsres@gsLineColor      = "black"
           gsres@gsLineThicknessF = 1.5

           do ipr = 0,npr-1
           do ilev = 0,nyy-1

              lnxx = (/prfPoly(ipr,ilev,0),prfPoly(ipr,ilev,1)/)
              lnyy = (/yy(ilev),yy(ilev)/)

              if (.not.all(ismissing(lnxx))) then
                 plydum(ilev_opt,iregion,ipr,ilev) = gsn_add_polyline(wks,xyplot(ilev_opt,iregion),lnxx,lnyy,gsres)
              end if

           end do
           end do
           delete(prfPoly)
           ;----------------------------------

           delete(xx)
           delete(yy)
           delete(res@trYMaxF)
           delete(res@trYMinF)
           delete(res@xyExplicitLegendLabels)
           delete(res@xyMarkers)
           delete(res@xyMarkerColors)
           delete(res@xyLineColors)

     end do ; zooming in or not 
     end do ; source or vicinity regions

     resP = True
     resP@gsnPanelYWhiteSpacePercent = 5 
     resP@gsnPanelXWhiteSpacePercent = 1 

     gsn_panel(wks,ndtooned(xyplot(:,:)),(/2,nregion/),resP)
     delete(wks)
     delete(resP)

   end if ; skip plotting the second case alone in a separate file
  end do ; case loop for profile plots
