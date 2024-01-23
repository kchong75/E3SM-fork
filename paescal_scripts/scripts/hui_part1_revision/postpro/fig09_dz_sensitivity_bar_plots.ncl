load "./utils/print_utils.ncl"

  ;-------------------------------------
  ; constants used for unit conversion
  ;-------------------------------------
  pi             = 3.14159265358979323846  ; from E3SM
  earth_radius   = 6.37122e6               ; from E3SM
  earth_sfc_area = 4.*pi*earth_radius^2

  kg2Tg  = 1e-9
  sec2yr = 86400.*365
  year2day = 365

  FillValue = -999.
  ;-------------------------------------

  itime = 0
  yearS = 2000
  yearE = 2009
  nyr   = yearE - yearS + 1

  dataFileTime = "ANN"

  ;----------------------------------------------------
  l_debug = True  ; print more information on screen
  l_debug = False ; print less information on screen

  dtime = 1800.

  dataFileRoot = "/pscratch/sd/h/huiwan/cflx/v1_papers/v1_runs/" 

  ;------------------------------------------
  ; case groups: different coupling schemes
  ;------------------------------------------
  cflx       = (/1, 2/)
  groupname  = (/"Original", "Revised"/)
  ngroup     = dimsizes(groupname)

  ;-----------------------------------------------
  ; cases in each group: different vertical grids
  ;-----------------------------------------------
  casenlev  = (/72,30/)
  ncase = dimsizes(casenlev)
  if (ncase.ne.2) then
     print("The script was written to compare 2 cases (simulations), but found ncase /= 2. Abort.")
  end if

  ;--------------------------------------------------------
  ; Aerosol species to process and their names in MAM
  ;--------------------------------------------------------
  species_longname = (/"Dust"/)
  species          = (/ "dst"/)
  nspecies         = dimsizes(species)

  ;-------------------------------------------------------------------------------------------------------------------
  ; E3SM's output variables used in this script are the interstitial mass burden and tendencies from CondiDiag.
  ;
  ; Variable names to be read from model output files have the pattern
  ;    QoI_prefix + QoI_name + chkpt
  ; where QoI_name = species + mode_id (e.g., "dst_a1", see mam4_species_and_modes.ncl)
  ;-------------------------------------------------------------------------------------------------------------------
  chkpt_label = (/"Sfc emis.", "Dry removal", "Burden"/)
  nchkpt      = dimsizes(chkpt_label)

  QoI_prefix = chkpt_label     ; create an array of the correct shape and type
  QoI_prefix = "cnd01_"        ; CondiDiag budget diagnostics

  ;-------------------------------------------------------------------------------
  ; The different regions to calculate average for, distinguished by emission.
  ; Grid columns chosen in each region are those with
  ;     emis_min < emis <=  emis_max
  ; Note that emis_max = 10. basically means no upper bound;
  ;           emis_min = -1. means including cells with no emission. 
  ;-------------------------------------------------------------------------------
  regionlbl = (/"Source-region total","Global total"/)
  emis_min = (/ 0.,  -1./)
  emis_max = (/10.,  10./)
  nregion  = dimsizes(emis_max)

  ;--------------------------------------------------------------------------------------
  ; Array for saving the spatial averages of all species, processes, regions, and years.
  ; At the end of the script execution, the multi-year mean and stddev values are
  ; presented in different combinations in different tables.
  ;--------------------------------------------------------------------------------------
  all_regionavg = new( (/ngroup,nspecies,nregion,nchkpt,ncase,nyr/), "float" )
      regionavg = new(        (/nspecies,nregion,nchkpt,ncase,nyr/), "float" )

  ;-------------------------------------------
  do igroup = 0,ngroup-1

     casename   = "cflx_"+cflx(igroup)+"_L"+casenlev+"_F20TRC5-CMIP6"

     emis_chkpt = casename ; make it an array
     emis_chkpt = "_v_CFLX"+cflx(igroup)+"_inc"

     print("")
     print("==============================")
     print(" simulations being compared:")
     print(" "+casename)
     print("==============================")

     ;------------------------ 
     ; checkpoints 

     chkpt = (/(/emis_chkpt(0), "_v_AERDRYRM_inc", "_v_AERWETRM"/) \
              ,(/emis_chkpt(1), "_v_AERDRYRM_inc", "_v_AERWETRM"/) \
              /)

     ; Unit conversion factors.
     ;  -  1/dtime is for converting increments to tendencies.
     ;  -  1/sec2yr is for canceling out a factor of sec2yr multiplied later.

     inc_fac =  (/ 1./dtime, 1./dtime, 1./sec2yr/)


     ;#############################################################################################
     ; Calculate regional averages for each species, year, QoI, and case
     ;#############################################################################################
     ;----------------------------------------
     ; Open data files and get some grid info
     ;----------------------------------------
     dataFilePath = dataFileRoot + casename + "/climo/"
     dataFileName = casename+".cam.h0."+yearS+"."+dataFileTime+".nc"
     dataFile = addfiles(dataFilePath+dataFileName,"r")

     area = tofloat(dataFile[0]->area)  ; cell area to be used in weighted averaging or global integral
     ncol = dimsizes(area)              ; size of horizonal dimension

     ;---------------------------------------------------
     ; Start reading and processing data
     ;---------------------------------------------------
     print("")
     do isp = 0,nspecies-1

        print("")
        print("============================")
        print(" SPECIES = "+species(isp))
        print("============================")
        print("")

        ; Get the list of tracer names for this species. Interstitial aerosols only.
        ; List is saved in variable QoI_name.

        load "./utils/mam4_species_and_modes.ncl"
        nqoi = dimsizes(QoI_name)  

        do iyr = 0,nyr-1

           year  = iyr + yearS
           dataFileName = casename+".cam.h0."+year+"."+dataFileTime+".nc"
           dataFile = addfiles(dataFilePath+dataFileName,"r")
           print("dataFile = "+dataFileName)

           do ichkpt = 0,nchkpt-1
 
              ;---------------------------
              ; get the sum over all modes 
              ;---------------------------
              ; initialize sum of all modes

              var2d = new( (/ncol,ncase/),"float")
              var2d = 0. 

              do iqoi = 0,nqoi-1
     
                 do icase = 0,ncase-1
                    varname = QoI_prefix(ichkpt)+QoI_name(iqoi)+chkpt(icase,ichkpt)
                    if (isfilevar(dataFile[icase],varname)) then
                      var2d(:,icase) = var2d(:,icase) + dataFile[icase]->$varname$(itime,:)*inc_fac(ichkpt)
                    end if
                 end do

                ;------------------ 
                if (l_debug) then 
                   if (isfilevar(dataFile[0],varname)) then
                      print("Done reading "+varname)
                   else
                      print(varname+" not found. Skip reading.")
                   end if
                end if
                ;------------------ 

              end do

              ;--------------------------------
              ; Save emissions for both cases
              ;--------------------------------
              if (chkpt_label(ichkpt).eq."Sfc emis.") then
                 emissions_global = var2d  ; both cases saved separately

                 if (l_debug) then
                    print("emissions saved. Min, max = "+min(emissions_global)+", "+max(emissions_global))
                 end if
              end if

              ;---------------------------------------------------------------------------------------
              ; Use emissions to select grid columns, then calculate regional averages/integrals.
              ; Note that selection is based on emissions in the specific case (simulation) and year.
              ;---------------------------------------------------------------------------------------
              do icase = 0,ncase-1
              do iregion = 0,nregion-1

                 idx = ind( emissions_global(:,icase).gt.emis_min(iregion) .and. \
                            emissions_global(:,icase).le.emis_max(iregion) )

                 if (.not.(all(ismissing(idx)))) then
                    ;regionavg(isp,iregion,ichkpt,icase,iyr) = sum(var2d(idx,icase)*area(idx))/sum(area(idx))                                   ; average flux or burden per m2 
                     regionavg(isp,iregion,ichkpt,icase,iyr) = sum(var2d(idx,icase)*area(idx))/sum(area(:)) * earth_sfc_area * kg2Tg * sec2yr   ; integral over region
                 end if
                 delete(idx)
              end do
              end do
 
              delete(var2d) 
           end do ; checkpoint loop

        end do ; year loop

        delete(QoI_name)
     end do ;species

     all_regionavg(igroup,:,:,:,:,:) = regionavg 

     print("")
     print("")
     print("")
     print("")
     ;#############################################################################################
     ; Done calculating regional averages for each species, year, QoI, and case
     ;#############################################################################################

  end do ; group

  ;********************************************************************************************************
  ; Create a data array, copy values already calculated above and diagnose some additional numbers.
  ;********************************************************************************************************
  ; dimension note: regionavg(isp,iregion,ichkpt,icase,iyr)

  row_titles = (/"(a) Dry removal in source regions (Tg/yr)", \
                 "(b) Global burden (Tg)", \
                 "(c) Global mean lifetime (day)"/)

  nrow = dimsizes(row_titles)

  isp = 0

  data_to_present = new( (/nrow,ngroup,ncase,nyr/), "float" )

  ii = -1

  ; dry removal

  ii = ii+1
  ichkpt = ind(chkpt_label.eq."Dry removal")
  iregion = 0  ; source region
  data_to_present(ii,:,:,:) = all_regionavg(:,isp,iregion,ichkpt,:,:)

  ; burden 

  ii = ii+1
  i9 = ii
  ichkpt = ind(chkpt_label.eq."Burden")
  iregion = 1  ; global 
  data_to_present(ii,:,:,:) = all_regionavg(:,isp,iregion,ichkpt,:,:)

  ; lifetime

  ichkpt_emis = ind(chkpt_label.eq."Sfc emis.")
  iregion = 1  ; global 
  total_source_tmp = all_regionavg(:,isp,iregion,ichkpt_emis,:,:)

  ii = ii+1
  data_to_present(ii,:,:,:) = data_to_present(i9,:,:,:)/total_source_tmp * year2day


  ;********************************************************************************************************
  ; Plotting
  ;********************************************************************************************************
  ; The plot file

  plotFileTimeTag = yearS+"-"+yearE+"_"+dataFileTime
  plotPath = "/global/cfs/projectdirs/m4359/www/huiwan/cflx/2023_rerun_for_gmd/"
  plotFile = plotPath+species(isp)+"_bar_plot_"+plotFileTimeTag
  plotFmt  = "pdf"

  wks = gsn_open_wks(plotFmt,plotFile)
  system("chmod a+r "+plotFile+"."+plotFmt)
  gsn_define_colormap(wks,"batlow")
  gsn_reverse_colormap(wks)

  print("plot file is "+plotFile+"."+plotFmt)
  if (.not.l_debug) then
     print("If you also want to have the plotted numbers printed in LaTex table format, set l_debug = True") 
  end if

  res               = True

  res@tmYLOn        = True
  res@tmYROn        = False
  res@tmYLBorderOn  = True
  res@tmYRBorderOn  = True

  res@tmXTOn        = False
  res@tmXTBorderOn  = True
  res@tmXBOn        = False
  res@tmXBBorderOn  = True

  res@tiXAxisString = " "
  res@tiYAxisString = " "

  FontHeightF = 0.022
  res@tiMainFontHeightF          = FontHeightF *1.1
  res@gsnCenterStringFontHeightF = FontHeightF
  res@tiXAxisFontHeightF         = FontHeightF *0.8
  res@tiYAxisFontHeightF         = FontHeightF *0.8
  res@tmXBLabelFontHeightF       = FontHeightF *0.8

  res@vpWidthF      = 0.6  ; change aspect ratio
  res@vpHeightF     = 0.6 

  res@gsnDraw       = False ; don't draw or advance the frame
  res@gsnFrame      = False

  res@tiMainFont      = "Helvetica-bold" 
  res@gsnCenterString = " "

  ; polygon resources

  pgres                  = True
  pgres@gsLineThicknessF = 1.5   ; make lines thicker
  pgres@gsLineDashPattern = 0

  ; line resources

  lnres = True
  lnres@gsLineThicknessF  = 2
  lnres@gsLineDashPattern = 0
  lnres@gsLineColor = 250 

  ; text resources

  txtres = True
  txtres@txFontHeightF = FontHeightF *0.9
  ;------

  plot           = new( nrow,"graphic" )
  dumMeanOutline = new( (/nrow,ngroup,ncase/),"graphic")
  dumMeanPolygon = new( (/nrow,ngroup,ncase/),"graphic")
  dumStdRange    = new( (/nrow,ngroup,ncase/),"graphic")
  dumStdLow      = new( (/nrow,ngroup,ncase/),"graphic")
  dumStdHgh      = new( (/nrow,ngroup,ncase/),"graphic")
  dumCaseLabel   = new( (/nrow,ngroup,ncase/),"graphic")
  dumGroupLabel  = new( (/nrow,ngroup      /),"graphic")

  ymax = (/ 3000.,40.,3./)*1.2
  ymin = (/    0.,0., 0./)

  do irow = 0,nrow-1  ; each "row" is a panel in figure

     res@tiMainString = row_titles(irow) 
     res@trYMaxF = ymax(irow) 
     res@trYMinF = ymin(irow) 

     res@trXMinF       =  0.
     res@trXMaxF       =  5.

     xleft   = 1.1
     dxGroup = 2.0
     dxCase  = 0.7
     dxBar   = 0.25
     dx      = dxBar

     ; Draw the frame first. The polygons and lines are added below in the nested loops.

     plot(irow) = gsn_csm_xy(wks,(/0,0/),(/res@trYMinF,res@trYMaxF/),res)

     do igroup = 0,ngroup-1
     do icase = 0,ncase-1

           zmean = abs(dim_avg   (data_to_present(irow,igroup,icase,:)))
           zstd  = dim_stddev(data_to_present(irow,igroup,icase,:)) 

           ; Show the mean value with a rectangle (polygon)

           x0 = xleft + igroup * dxGroup + icase * dxCase
           xp = (/x0-dx,x0+dx,x0+dx,x0-dx,x0-dx/)
           yp = (/0.,   0.,   zmean,zmean, 0. /)

           pgres@gsFillColor = (irow+1)*45
           pgres@gsFillIndex = (igroup)*3
           pgres@gsFillScaleF = 0.4

           dumMeanPolygon(irow,igroup,icase) = gsn_add_polygon (wks,plot(irow),xp,yp,pgres)
           dumMeanOutline(irow,igroup,icase) = gsn_add_polyline(wks,plot(irow),xp,yp,pgres)

           ; Show the 2-std range as whiskers.

           ; The range

           xln = (/x0,x0/)
           yln = (/zmean-zstd,zmean+zstd/)
           dumStdRange(irow,igroup,icase) = gsn_add_polyline(wks,plot(irow),xln,yln,lnres)

           ; Lines showing the bounds of the range

           xln = (/x0-0.5*dx,x0+0.5*dx/)

           yln = (/zmean-zstd,zmean-zstd/)
           dumStdLow(irow,igroup,icase) = gsn_add_polyline(wks,plot(irow),xln,yln,lnres)

           yln = (/zmean+zstd,zmean+zstd/)
           dumStdHgh(irow,igroup,icase) = gsn_add_polyline(wks,plot(irow),xln,yln,lnres)

           ; Add group label

           xt = x0 
           yt = res@trYMaxF *0.85
           txt = "L"+casenlev(icase)

           txtres@txFont = "Helvetica" 
           dumCaseLabel(irow,igroup,icase) = gsn_add_text(wks,plot(irow),txt,xt,yt,txtres)

        end do

        ; Add group label

        xt = xleft + igroup * dxGroup + 0.5*dxCase
        yt = res@trYMaxF *0.92
        txt = groupname(igroup)

        txtres@txFont = "Helvetica-bold" 
        dumGroupLabel(irow,igroup) = gsn_add_text(wks,plot(irow),txt,xt,yt,txtres)

     end do

  end do

  ;======================
  ; panel plots together
  ;======================
  pres = True
  pres@gsnPanelYWhiteSpacePercent = 5

  gsn_panel(wks,plot(:),(/1,nrow/),pres)


;@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
; Done with plotting. The rest of the script prints the plotted
; numbers as a LaTex table for verification.
;@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
if (l_debug) then

  do igroup = 0,ngroup-1

  ;--------------------------------------
  ; Multi-year mean and stddev
  ;--------------------------------------
  budget_terms_mean = dim_avg   (data_to_present(:,igroup,:,:))   ; ending dims: nrow,ncase
  budget_terms_std  = dim_stddev(data_to_present(:,igroup,:,:))   ; ending dims: nrow,ncase

  ; Std. dev. expressed as percentage of mean

  denom = where( budget_terms_mean.ne.0, budget_terms_mean, FillValue)
  denom@_FillValue = FillValue
  budget_terms_std_pct =  budget_terms_std/abs(denom) *100.
  budget_terms_std_pct = where( ismissing(budget_terms_std_pct), 0., budget_terms_std_pct )  ; change missing values to 0 for printing
  delete(denom)

  ;--------------------------------------
  ; Relative diff of multi-year mean
  ;--------------------------------------
  denom = where( budget_terms_mean(:,0).ne.0, budget_terms_mean(:,0), FillValue)
  denom@_FillValue = FillValue
  reldiff_of_mean = budget_terms_mean(:,1)/denom - 1.
  reldiff_of_mean = where( ismissing(reldiff_of_mean), 0., reldiff_of_mean )  ; change missing values to 0 for printing
  delete(denom)

  ;******************************************************
  ; Print table(s) of budget by species
  ;******************************************************
  species_table = new( nrow+2, "string")  ; the 2 extra rows are title rows
  species_table(:) = ""

  isp = 0
     explbl4table = groupname(igroup)+" L"+casenlev

     print("")
     print("=====================================================================")
     print(" Table for "+species_longname(isp)+", case group = "+groupname(igroup))
     print("=====================================================================")
     print("")

     print("\middlehline")

     ; first column of table

     species_table(0)  = "\multirow{ 2}{*}{\bf "+species_longname(isp)+" budget} &"
     species_table(1)  = " & "
     species_table(2:) = row_titles + " & "

     ; header row 0: different regions

     irwh = 0
     species_table(irwh) = species_table(irwh)+ " \multicolumn{3}{c}{{\bf region of choice }} &&"

     ; header row 1: 3 columns of numbers in each region

     irwh = irwh+1
     species_table(irwh) = species_table(irwh)+ " "+explbl4table(0)+" & "+explbl4table(1)+" & Rel. diff. &&"

     ; numbers

     irwc0 = irwh+1    ; first row that contains numbers

     do ii = 0,nrow-1

        if (ii.lt.nrow-2) then
           floatFmt = "%-.0f"     ; Process rates are large. Show them as integers. 
        else
           floatFmt = "%-.1f"     ; Burden and lifetime are small numbers. Show one digit after decimal point.
        end if

        pcntgFmt = "%-.0f"        ; Show percentages as integers

        str0a = sprintf_PlsMns_LaTex( floatFmt, budget_terms_mean(ii,0) )
        str1a = sprintf_PlsMns_LaTex( floatFmt, budget_terms_mean(ii,1) )
        str2  = sprintf_PlsMns_LaTex( pcntgFmt, reldiff_of_mean(ii)*100 )

        str0b = sprintf( pcntgFmt, budget_terms_std_pct(ii,0) )
        str1b = sprintf( pcntgFmt, budget_terms_std_pct(ii,1) )

        species_table(ii+irwc0) = species_table(ii+irwc0) + str0a+" $(\pm"+str0b+"\%)$" \
                                                    +" & "+ str1a+" $(\pm"+str1b+"\%)$" \
                                                    +" & " +str2+"\% &&"
     end do
     species_table(:) = species_table(:) + "\\"

     print(""+species_table)
     print("")

  end do  ; group

end if ; l_debug
;@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
