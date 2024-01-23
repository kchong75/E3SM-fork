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
  yearE = 2000
  nyr   = yearE - yearS + 1


  ;----------------------------------------------------
  l_debug = True  ; print more information on screen
  l_debug = False ; print less information on screen


  casename = (/"cflx_cpl_opt1_F20TRC5-CMIP6_pm-cpu","cflx_cpl_opt2_F20TRC5-CMIP6_pm-cpu"/)   ; revised vs original, L72
  explbl4table = (/"Original","Revised"/)

 ;casename = (/"compy_F20TRC5-CMIP6_ne30_ndg_cflx_cpl_opt1","compy_F20TRC5-CMIP6_ne30_ndg_cflx_cpl_opt1_L30"/)  ; L72 vs L30, original
 ;explbl4table = (/"L72","L30"/)

 ;casename = (/"compy_F20TRC5-CMIP6_ne30_ndg_cflx_cpl_opt2","compy_F20TRC5-CMIP6_ne30_ndg_cflx_cpl_opt2_L30"/)  ; L72 vs L30, revised
 ;explbl4table = (/"L72","L30"/)

  ncase = dimsizes(casename)
  if (ncase.ne.2) then
     print("The script was written to compare 2 cases (simulations), but found ncase /= 2. Abort.")
  end if
  
  dataFilePath = "/pscratch/sd/h/huiwan/cflx/v1_papers/v1_runs/" + casename + "/climo/"
  ;------------------------
  print("")
  print("=====================")
  print(" simulations being compared:")
  print(" "+casename)
  print("=====================")

  ;--------------------------------------------------------
  ; Aerosol species to process and their names in MAM
  ;--------------------------------------------------------
  species_longname = (/"Dust","Sea salt","MOA","BC","POA"/)
  species          = (/ "dst",     "ncl","mom","bc","pom"/)
  nspecies         = dimsizes(species)

  ;-----------------------------------------------------------------------------------------------------
  ; The model output variables used in this script are
  ;   - aerosol budget diagnostics from MAM (as these are already available even without CondiDiag), and
  ;   - interstitial mass burden from CondiDiag.
  ;
  ; Variable names to be read from model output file have the pattern
  ;    QoI_prefix + QoI_name + chkpt
  ; where QoI_name = species + mode_id (e.g., "dst_a1", see mam4_species_and_modes.ncl)
  ;------------------------------------------------------------------------------------
  QoI_prefix = (/"SF",                  "",           "",           "", "cnd01_"/)
  chkpt      = (/"",                 "DDF",   "_mixnuc1",      "SFWET", "_v_AERWETRM"/)
  inc_fac    = (/ 1.,                  -1.,           1.,           1., 1./sec2yr/)    ; for sign and unit conversions
  label      = (/"Sfc emis.","Dry removal","Act./resus.","Wet removal", "Burden"/)     ; text in LaTex tables

  nchkpt = dimsizes(chkpt)

  ;-------------------------------------------------------------------------------
  ; The different regions to calculate average for, sorted by emission.
  ; Grid cells chosen in each region are those with
  ;     emis_min < emis <=  emis_max
  ; Note that emis_max = 10. basically means no upper bound;
  ;           emis_min = -1. means including cells with no emission. 
  ;-------------------------------------------------------------------------------
  regionlbl = (/"Source-region mean","Non-source-region mean","Global mean"/)
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

  area = tofloat(dataFile[0]->area)  ; cell area to be used in weighted averaging
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

     load "./utils/mam4_species_and_modes.ncl"
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
                print("Done reading "+varname)
             else
                print(varname+" not found. Skip reading.")
             end if
  
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

           ;--------------------------------------------------------------
           ; Use emissions to sort grid columns, then calculate averages 
           ;--------------------------------------------------------------
           do icase = 0,ncase-1
           do iregion = 0,nregion-1

              idx = ind( emissions_global(:,icase).gt.emis_min(iregion) .and. \
                         emissions_global(:,icase).le.emis_max(iregion) )

              if (.not.(all(ismissing(idx)))) then
                  regionavg(isp,iregion,ichkpt,icase,iyr) = sum(var2d(idx,icase)*area(idx))/sum(area(idx)) * earth_sfc_area * kg2Tg * sec2yr
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

        iregion = 2 ; global
        do ichkpt = 0,nchkpt-1
        print("Global mean: "+sprintf("%10.2f",regionavg_mean(isp,iregion,ichkpt,0))+" Tg/yr vs "\
                             +sprintf("%10.2f",regionavg_mean(isp,iregion,ichkpt,1))+" Tg/yr from "+label(ichkpt) )
        end do

        do ichkpt = 0,nchkpt-1
        print("Global std : "+sprintf("%10.2f",regionavg_std(isp,iregion,ichkpt,0))+" Tg/yr vs "\
                             +sprintf("%10.2f",regionavg_std(isp,iregion,ichkpt,1))+" Tg/yr from "+label(ichkpt) )
        end do

        delete(regionavg_mean)
        delete(regionavg_std)
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

  row_titles = (/"Emissions (Tg yr$^{-1}$)", \
                 "Dry removal (Tg yr$^{-1}$)", \
                 "Activation (Tg yr$^{-1}$)",  \
                 "Wet removal (Tg yr$^{-1}$)", \
                 "Resolved transport (Tg yr$^{-1}$)", \
                 "Burden (Tg)", \
                 "Lifetime (day)"/)

  budget_terms = new( (/nspecies,nregion,dimsizes(row_titles),ncase,nyr/), "float" )
  ii = -1

  ; emissions

  ii = ii+1
  ichkpt = ind(QoI_prefix.eq."SF")
  budget_terms(:,:,ii,:,:) = regionavg(:,:,ichkpt,:,:)

  ; dry removal

  ii = ii+1
  ichkpt = ind(chkpt.eq."DDF")
  budget_terms(:,:,ii,:,:) = regionavg(:,:,ichkpt,:,:)

  ; activation

  ii = ii+1
  ichkpt = ind(chkpt.eq."_mixnuc1")
  budget_terms(:,:,ii,:,:) = ( regionavg(:,:,ichkpt,:,:) )

  ; wet removal total: in-cloud + below-cloud

  ii = ii+1
  ichkpt = ind(chkpt.eq."SFWET")
  budget_terms(:,:,ii,:,:) = ( regionavg(:,:,ichkpt,:,:) )

  ; abs(transport) = (approximately) surface emission - abs( dry + wet + act )

  ii = ii+1
  i0 = ii
  ichkpt_emis = ind(QoI_prefix.eq."SF")
  budget_terms(:,:,ii,:,:) = -   regionavg(:,:,ichkpt_emis,             :,:) \
                             - ( regionavg(:,:,ind(chkpt.eq."DDF"),     :,:) \
                                +regionavg(:,:,ind(chkpt.eq."_mixnuc1"),:,:) \
                                +regionavg(:,:,ind(chkpt.eq."SFWET"),   :,:) )

  ; burden 

  ii = ii+1
  i9 = ii
  ichkpt = ind(chkpt.eq."_v_AERWETRM")
  budget_terms(:,:,i9,:,:) = regionavg(:,:,ichkpt,:,:)

  ; lifetime
  ;                          glb or src region mean,               emissions as src,                    transport as src
  total_source_tmp = where(  regionavg(:,:,ichkpt_emis,:,:).gt.0., regionavg(:,:,ichkpt_emis,:,:), budget_terms(:,:,i0,:,:) )

  ii = ii+1
  budget_terms(:,:,ii,:,:) = budget_terms(:,:,i9,:,:)/total_source_tmp *365.
 
  ;--------------------------------------
  ; Multi-year mean and stddev
  ;--------------------------------------
  budget_terms_mean = dim_avg   (budget_terms)
  budget_terms_std  = dim_stddev(budget_terms)

  ;--------------------------------------
  ; Relative diff of multi-year mean
  ;--------------------------------------
  denom = where( budget_terms_mean(:,:,:,0).ne.0, budget_terms_mean(:,:,:,0), FillValue)
  denom@_FillValue = FillValue
  reldiff_of_mean = budget_terms_mean(:,:,:,1)/denom - 1.
  reldiff_of_mean = where( ismissing(reldiff_of_mean), 0., reldiff_of_mean )
  delete(denom)

  ;******************************************************
  ; Print table(s) to budget(s) by species
  ;******************************************************
  species_table = new( dimsizes(row_titles)+2, "string")
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
     do ii = 0,dimsizes(row_titles)-1
     do iregion = 0,nregion-1

        floatFmt = "%-.1f"
        pcntgFmt = "%-.0f"
        str0 = sprintf_PlsMns_LaTex( floatFmt, budget_terms_mean(isp,iregion,ii,0) )
        str1 = sprintf_PlsMns_LaTex( floatFmt, budget_terms_mean(isp,iregion,ii,1) )
        str2 = sprintf_PlsMns_LaTex( pcntgFmt, reldiff_of_mean(isp,iregion,ii)*100 )

        species_table(ii+irw) = species_table(ii+irw) + str0+" & "+ str1+" & " +str2+"\% &"
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
     iprocE = dimsizes(row_titles)-1
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

        floatFmt = "%-.1f"
        pcntgFmt = "%-.0f"
        str0 = sprintf_PlsMns_LaTex( floatFmt, budget_terms_mean(isp,iregion,ii,0) )
        str1 = sprintf_PlsMns_LaTex( floatFmt, budget_terms_mean(isp,iregion,ii,1) )
        str2 = sprintf_PlsMns_LaTex( pcntgFmt, reldiff_of_mean(isp,iregion,ii)*100 )

        table_row(isp+irw) = table_row(isp+irw) + str0+" & "+ str1+" & " +str2+"\% &"

     end do
     end do
     table_row(:) = table_row(:) + "\\"

     print(""+table_row)
     print("")
     delete(table_row)
  end do

