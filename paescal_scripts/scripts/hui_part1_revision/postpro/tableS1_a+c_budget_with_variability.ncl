load "./utils/print_utils.ncl"

  ;-------------------------------------
  ; constants used for unit conversion
  ;-------------------------------------
  pi             = 3.14159265358979323846  ; from E3SM
  earth_radius   = 6.37122e6               ; from E3SM
  earth_sfc_area = 4.*pi*earth_radius^2

  kg2Tg  = 1e-9
  sec2yr = 86400.*365

  FillValue = -999.
  ;-------------------------------------

  itime = 0
  yearS = 2000
  yearE = 2009
  nyr   = yearE - yearS + 1

  ;----------------------------------------------------
  l_debug = True  ; print more information on screen
  l_debug = False ; print less information on screen


  casename = (/"cflx_1_L72_F20TRC5-CMIP6","cflx_2_L72_F20TRC5-CMIP6"/)   ; revised vs original, L72
  explbl4table = (/"Original","Revised"/)

  ncase = dimsizes(casename)
  if (ncase.ne.2) then
     print("The script was written to compare 2 cases (simulations), but found ncase /= 2. Abort.")
  end if
  
  dataFilePath = "/pscratch/sd/h/huiwan/cflx/v1_papers/v1_runs/" + casename + "/climo/"
  ;------------------------
  print("")
  print("====================================")
  print(" simulations being compared:")
  print(" "+casename)
  print("====================================")

  ;--------------------------------------------------------
  ; Aerosol species to process and their names in MAM
  ;--------------------------------------------------------
  species_longname = (/"Dust","Sea salt","MOA","BC","POA"/)
  species          = (/ "dst",     "ncl","mom","bc","pom"/)
  nspecies         = dimsizes(species)

  ;-----------------------------------------------------------------------------------------------------
  ; E3SM's output variables used in this script are
  ;   - aerosol budget diagnostics from MAM, and
  ;   - interstitial mass burden from CondiDiag.
  ;
  ; Variable names to be read from model output files have the pattern
  ;    QoI_prefix + QoI_name + chkpt
  ; where QoI_name = species + mode_id (e.g., "dst_a1", see mam4_species_and_modes.ncl)
  ;------------------------------------------------------------------------------------
  QoI_prefix = (/"SF",                  "",           ""/)
  chkpt      = (/"",                 "DDF",      "SFWET"/)
  inc_fac    = (/ 1.,                  -1.,           1./)    ; for sign and unit conversions
  label      = (/"Sfc emis.","Dry removal","Wet removal"/)    ; text printed on screen

  nchkpt = dimsizes(chkpt)

  ;-------------------------------------------------------------------------------
  ; The different regions to calculate average for, sorted by emission.
  ; Grid cells chosen in each region are those with
  ;     emis_min < emis <=  emis_max
  ; Note that emis_max = 10. basically means no upper bound;
  ;           emis_min = -1. means including cells with no emission. 
  ;-------------------------------------------------------------------------------
  regionlbl = (/"Source-region total","Non-source-region total","Global total"/)
  emis_min = (/ 0., -1., -1./)
  emis_max = (/10.,  0., 10./)
  nregion  = dimsizes(emis_max)

  ;--------------------------------------------------------------------------------------
  ; Array for saving the spatial averages of all species, processes, regions, and years.
  ; At the end of the script execution, the multi-year mean and stddev values are
  ; presented in different combinations in different tables.
  ;--------------------------------------------------------------------------------------
  regionavg = new( (/nspecies,nregion,nchkpt,ncase,nyr/), "float" )

  ;----------------------------------------
  ; Open data file and get some grid info
  ;----------------------------------------
  dataFileName = casename+".cam.h0."+yearS+".ANN.nc"
  dataFile = addfiles(dataFilePath+dataFileName,"r")

  area = tofloat(dataFile[0]->area)  ; cell area to be used in weighted averaging or global integral
  ncol = dimsizes(area)              ; size of horizonal dimension

  ;---------------------------------------------------
  ; Start processing reading and processing data
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

     load "./utils/mam4_species_and_modes_ax_and_cx.ncl"  ; interstitial + cloudborne
     nqoi = dimsizes(QoI_name)  

     do iyr = 0,nyr-1

        year  = iyr + yearS
        dataFileName = casename+".cam.h0."+year+".ANN.nc"
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
  
             varname = QoI_prefix(ichkpt)+QoI_name(iqoi)+chkpt(ichkpt)
  
             if (isfilevar(dataFile[0],varname)) then
                do icase = 0,ncase-1
                   var2d(:,icase) = var2d(:,icase) + dataFile[icase]->$varname$(itime,:)*inc_fac(ichkpt)
                end do
             end if

            ;;------------------ 
            ;if (l_debug) then 
            ;   if (isfilevar(dataFile[0],varname)) then
            ;      print("Done reading "+varname)
            ;   else
            ;      print(varname+" not found. Skip reading.")
            ;   end if
            ;end if
            ;;------------------ 

           end do

           ;--------------------------------
           ; Save emissions for both cases
           ;--------------------------------
           if (QoI_prefix(ichkpt).eq."SF") then
              emissions_global = var2d  ; both cases saved separately

              if (l_debug) then
                 print("emissions saved. Min, max = "+min(emissions_global)+", "+max(emissions_global))
              end if
           end if

           ;---------------------------------------------------------------------------------------
           ; Use emissions to select grid columns, then calculate averages/integrals
           ; Note that selection is based on emissions in the specific case and year.
           ;---------------------------------------------------------------------------------------
           do icase = 0,ncase-1
           do iregion = 0,nregion-1

              idx = ind( emissions_global(:,icase).gt.emis_min(iregion) .and. \
                         emissions_global(:,icase).le.emis_max(iregion) )

              if (.not.(all(ismissing(idx)))) then
                 ;regionavg(isp,iregion,ichkpt,icase,iyr) = sum(var2d(idx,icase)*area(idx))/sum(area(idx)) * earth_sfc_area * kg2Tg * sec2yr
                  regionavg(isp,iregion,ichkpt,icase,iyr) = sum(var2d(idx,icase)*area(idx))/sum(area(:)) * earth_sfc_area * kg2Tg * sec2yr
              end if
              delete(idx)
           end do
           end do
 
           delete(var2d) 
        end do ; checkpoint loop

     end do ; year loop

     ;------------------------
     if (l_debug) then
        print("")
        regionavg_mean = dim_avg   (regionavg)
        regionavg_std  = dim_stddev(regionavg)

        FillValue = -999.
        denom = where( abs(regionavg_mean).gt.1e-36, regionavg_mean, FillValue )
        denom@_FillValue = FillValue
        regionavg_std_over_mean  = regionavg_std/abs(denom)*100.
        delete(denom)

        iregion = 2 ; global
        do ichkpt = 0,nchkpt-1
        print("Global mean: "+sprintf("%10.2f",regionavg_mean(isp,iregion,ichkpt,0))+" Tg/yr vs "\
                             +sprintf("%10.2f",regionavg_mean(isp,iregion,ichkpt,1))+" Tg/yr from "+label(ichkpt) )
        end do

        do ichkpt = 0,nchkpt-1
        print("Global std : "+sprintf("%10.2f",regionavg_std(isp,iregion,ichkpt,0))+" Tg/yr vs "\
                             +sprintf("%10.2f",regionavg_std(isp,iregion,ichkpt,1))+" Tg/yr from "+label(ichkpt) )
        end do

        print("")
        do ichkpt = 0,nchkpt-1
        print("Global std/mean : "+sprintf("%10.2f",regionavg_std_over_mean(isp,iregion,ichkpt,0))+" % "\
                                  +sprintf("%10.2f",regionavg_std_over_mean(isp,iregion,ichkpt,1))+" % from "+label(ichkpt) )
        end do

        print("")
        iregion = 0 ; source region
        do ichkpt = 0,nchkpt-1
        print("Src region std/mean : "+sprintf("%10.2f",regionavg_std_over_mean(isp,iregion,ichkpt,0))+" % "\
                                      +sprintf("%10.2f",regionavg_std_over_mean(isp,iregion,ichkpt,1))+" % from "+label(ichkpt) )
        end do

        delete(regionavg_mean)
        delete(regionavg_std)
       ;exit
     end if
     ;------------------------

     delete(QoI_name)
  end do ;species

  print("")
  print("")
  print("")
  print("")
  ;********************************************************************************************************
  ; Create a data array, copy values already calculated above and diagnose some additional numbers.
  ;********************************************************************************************************
  ; dimension note: regionavg(isp,iregion,ichkpt,icase,year)

  row_titles = (/"Emissions",  \
                 "Dry removal",\
                 "Wet removal" /)

  nrow = dimsizes(row_titles)

  budget_terms = new( (/nspecies,nregion,nrow,ncase,nyr/), "float" )
  ii = -1

  ; emissions

  ii = ii+1
  ichkpt = ind(QoI_prefix.eq."SF")
  budget_terms(:,:,ii,:,:) = regionavg(:,:,ichkpt,:,:)

  ; dry removal

  ii = ii+1
  ichkpt = ind(chkpt.eq."DDF")
  budget_terms(:,:,ii,:,:) = regionavg(:,:,ichkpt,:,:)

  ; wet removal total: in-cloud + below-cloud

  ii = ii+1
  ichkpt = ind(chkpt.eq."SFWET")
  budget_terms(:,:,ii,:,:) = ( regionavg(:,:,ichkpt,:,:) )

  ;--------------------------------------
  ; Multi-year mean and stddev
  ;--------------------------------------
  budget_terms_mean = dim_avg   (budget_terms)
  budget_terms_std  = dim_stddev(budget_terms)

  ; Std. dev. expressed as percentage of mean

  denom = where( budget_terms_mean.ne.0, budget_terms_mean, FillValue)
  denom@_FillValue = FillValue
  budget_terms_std_pct  =  budget_terms_std/abs(denom) *100.
  budget_terms_std_pct = where( ismissing(budget_terms_std_pct), 0., budget_terms_std_pct )
  delete(denom)

  ;--------------------------------------
  ; Relative diff of multi-year mean
  ;--------------------------------------
  denom = where( budget_terms_mean(:,:,:,0).ne.0, budget_terms_mean(:,:,:,0), FillValue)
  denom@_FillValue = FillValue
  reldiff_of_mean = budget_terms_mean(:,:,:,1)/denom - 1.
  reldiff_of_mean = where( ismissing(reldiff_of_mean), 0., reldiff_of_mean )
  delete(denom)

  ;******************************************************
  ; Print table(s) of budget by species
  ;******************************************************
  species_table = new( nrow+2, "string")  ; the 2 extra rows are title rows
  species_table(:) = ""

  if (.not.l_debug) then
     nspectmp = 0
  else
     nspectmp = dimsizes(species)-1
  end if

  do isp = 0,nspectmp

     print("")
     print("=========================================")
     print(" Budget table for "+species_longname(isp))
     print("=========================================")
     print("")

     print("\middlehline")
     species_table(0)  = "\multirow{ 2}{*}{\bf "+species_longname(isp)+" budget} &"
     species_table(1)  = " & "
     species_table(2:) = row_titles + " & "

     ; title row 0

     irw = 0
     do iregion = 0,nregion-1
        species_table(irw) = species_table(irw)+ " \multicolumn{3}{c}{{\bf "+regionlbl(iregion)+"}} &"
     end do

     irw = irw+1
     do iregion = 0,nregion-1
        species_table(irw) = species_table(irw)+ " "+explbl4table(0)+" & "+explbl4table(1)+" & Rel. diff. &"
     end do

     ; results

     irw = irw+1
     do ii = 0,nrow-1
     do iregion = 0,nregion-1

        floatFmt = "%-.2f"
        pcntgFmt = "%-.0f"

        str0a = sprintf_PlsMns_LaTex( floatFmt, budget_terms_mean(isp,iregion,ii,0) )
        str1a = sprintf_PlsMns_LaTex( floatFmt, budget_terms_mean(isp,iregion,ii,1) )
        str2  = sprintf_PlsMns_LaTex( pcntgFmt, reldiff_of_mean(isp,iregion,ii)*100 )

        str0b = sprintf( pcntgFmt, budget_terms_std_pct(isp,iregion,ii,0) )
        str1b = sprintf( pcntgFmt, budget_terms_std_pct(isp,iregion,ii,1) )

        species_table(ii+irw) = species_table(ii+irw) + str0a+" $(\pm "+str0b+"\%)$" \
                                                +" & "+ str1a+" $(\pm "+str1b+"\%)$" \
                                                +" & " +str2+"\% &"
     end do
     end do
     species_table(:) = species_table(:) + "\\"

     print(""+species_table)
     print("")
  end do


  ;******************************************************
  ; Print table showing budgets by process
  ;******************************************************

  if (.not.l_debug) then
     iprocS = 1 ; dry removal
     iprocE = 1 ; dry removal
  else
     iprocS = 0
     iprocE = dimsizes(row_titles)-3
  end if

  do ii = iprocS,iprocE
     print("============================================================")
     print(" Process rate table of "+ row_titles(ii))
     print("============================================================")


     table_row = new( nspecies+2, "string")
     table_row(:) = ""
     print("\middlehline")

     table_row(0)  = "\multirow{ 2}{*}{\bf "+row_titles(ii)+"} & "
     table_row(1)  = "& "
     table_row(2:) = species_longname + " & "

     ; title row 0: lists regions
     irw = 0 
     do iregion = 0,nregion-1
        table_row(irw) = table_row(irw)+ " \multicolumn{3}{c}{\bf "+regionlbl(iregion)+"} &"
     end do

     ; title row 1: for each region, we have 2 cases and a diff.
     irw = irw+1
     do iregion = 0,nregion-1
        table_row(irw) = table_row(irw)+ " "+explbl4table(0)+" & "+explbl4table(1)+" & Rel. diff. &"
     end do

     ; results

     irw = irw+1

     do isp = 0,nspecies-1
     do iregion = 0,nregion-1

        if (any(species(isp).eq.(/"dst","ncl"/))) then
           floatFmt = "%-.0f"
        else
           floatFmt = "%-.2f"
        end if
        pcntgFmt = "%-.0f"

        str0a = sprintf_PlsMns_LaTex( floatFmt, budget_terms_mean(isp,iregion,ii,0) )
        str1a = sprintf_PlsMns_LaTex( floatFmt, budget_terms_mean(isp,iregion,ii,1) )
        str2  = sprintf_PlsMns_LaTex( pcntgFmt, reldiff_of_mean(isp,iregion,ii)*100 )

        str0b = sprintf( pcntgFmt, budget_terms_std_pct(isp,iregion,ii,0) )
        str1b = sprintf( pcntgFmt, budget_terms_std_pct(isp,iregion,ii,1) )

        table_row(isp+irw) = table_row(isp+irw) +       str0a+" $(\pm "+str0b+"\%)$" \
                                                +" & "+ str1a+" $(\pm "+str1b+"\%)$" \
                                                +" & " +str2+"\% &"

     end do
     end do
     table_row(:) = table_row(:) + "\\"

     print(""+table_row)
     print("")
     delete(table_row)
  end do

