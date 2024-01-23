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

  ;----------------------------------------------------
  l_debug = True  ; print more information on screen
  l_debug = False ; print less information on screen

  dtime = 1800.
  casename = (/"cflx_1_L30_F20TRC5-CMIP6","cflx_2_L30_F20TRC5-CMIP6"/)   ; revised vs original, L30
  casename = (/"cflx_1_L72_F20TRC5-CMIP6","cflx_2_L72_F20TRC5-CMIP6"/)   ; revised vs original, L72
  explbl4table = (/"Original","Revised"/)

  ncase = dimsizes(casename)
  if (ncase.ne.2) then
     print("The script was written to compare 2 cases (simulations), but found ncase /= 2. Abort.")
  end if
  
  dataFilePath = "/pscratch/sd/h/huiwan/cflx/v1_papers/v1_runs/" + casename + "/climo/"
  ;------------------------
  print("")
  print("============================")
  print(" simulations being compared:")
  print(" "+casename)
  print("============================")

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
  label      = (/"Sfc emis.",        "Dry removal",      "Act./resus.",  "Wet removal",     "Transport",       "Burden"/)

  QoI_prefix = label     ; create an array of the correct shape and type
  QoI_prefix = "cnd01_"  ; CondiDiag budget diagnostics

  chkpt      = (/(/"_v_CFLX1_inc",     "_v_AERDRYRM_inc",  "_v_STCLD_inc", "_v_AERWETRM_inc", "_v_PBCINI_inc", "_v_AERWETRM"/) \
                ,(/"_v_CFLX2_inc",     "_v_AERDRYRM_inc",  "_v_STCLD_inc", "_v_AERWETRM_inc", "_v_PBCINI_inc", "_v_AERWETRM"/) \
                /)

  ; Unit conversion factors.
  ;  -  1/dtime is for converting increments to tendencies.
  ;  -  1/sec2yr is for canceling out a factor of sec2yr multiplied later.

  inc_fac    =   (/ 1./dtime,           1./dtime,          1./dtime,        1./dtime,         1./dtime,        1./sec2yr/)

  nchkpt     = dimsizes(inc_fac)

  ;-------------------------------------------------------------------------------
  ; The different regions to calculate average for, distinguished by emission.
  ; Grid cells/columns chosen in each region are those with
  ;     emis_min < emis <=  emis_max
  ; Note that emis_max = 10. basically means no upper bound;
  ;           emis_min = -1. means including cells/columns with no emission. 
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
  ; Open data files and get some grid info
  ;----------------------------------------
  dataFileName = casename+".cam.h0."+yearS+".ANN.nc"
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
           if (label(ichkpt).eq."Sfc emis.") then
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

  print("")
  print("")
  print("")
  print("")

  ;********************************************************************************************************
  ; Create a data array, copy values already calculated above and diagnose some additional numbers.
  ;********************************************************************************************************
  ; dimension note: regionavg(isp,iregion,ichkpt,icase,iyr)

  row_titles = (/"Emissions [Tg yr$^{-1}$]", \
                 "Dry removal [Tg yr$^{-1}$]", \
                 "Activation [Tg yr$^{-1}$]",  \
                 "Wet removal [Tg yr$^{-1}$]", \
                 "Resolved transport [Tg yr$^{-1}$]", \
                 "Burden [Tg]", \
                 "Lifetime [day]"/)

  nrow = dimsizes(row_titles)

  budget_terms = new( (/nspecies,nregion,nrow,ncase,nyr/), "float" )
  ii = -1

  ; emissions

  ii = ii+1
  ichkpt = ind(label.eq."Sfc emis.")
  budget_terms(:,:,ii,:,:) = regionavg(:,:,ichkpt,:,:)
  ichkpt_emis = ichkpt   ; save for use later

  ; dry removal

  ii = ii+1
  ichkpt = ind(label.eq."Dry removal")
  budget_terms(:,:,ii,:,:) = regionavg(:,:,ichkpt,:,:)

  ; activation

  ii = ii+1
  ichkpt = ind(label.eq."Act./resus.")
  budget_terms(:,:,ii,:,:) = ( regionavg(:,:,ichkpt,:,:) )

  ; wet removal total: in-cloud + below-cloud

  ii = ii+1
  ichkpt = ind(label.eq."Wet removal")
  budget_terms(:,:,ii,:,:) = ( regionavg(:,:,ichkpt,:,:) )

  ; transport 

  ii = ii+1
  ichkpt = ind(label.eq."Transport")
  budget_terms(:,:,ii,:,:) = ( regionavg(:,:,ichkpt,:,:) )
  i0 = ii       ; save for use later

  ; burden 

  ii = ii+1
  ichkpt = ind(label.eq."Burden")
  budget_terms(:,:,ii,:,:) = regionavg(:,:,ichkpt,:,:)
  i9 = ii       ; save for use later

  ; lifetime
  ;                          glb or src region mean,               emissions as src,                    transport as src
  total_source_tmp = where(  regionavg(:,:,ichkpt_emis,:,:).gt.0., regionavg(:,:,ichkpt_emis,:,:), budget_terms(:,:,i0,:,:) )

  ii = ii+1
  budget_terms(:,:,ii,:,:) = budget_terms(:,:,i9,:,:)/total_source_tmp * year2day
 
  ;--------------------------------------
  ; Multi-year mean and stddev
  ;--------------------------------------
  budget_terms_mean = dim_avg   (budget_terms)
  budget_terms_std  = dim_stddev(budget_terms)

  ; Std. dev. expressed as percentage of mean

  denom = where( budget_terms_mean.ne.0, budget_terms_mean, FillValue)
  denom@_FillValue = FillValue
  budget_terms_std_pct =  budget_terms_std/abs(denom) *100.
  budget_terms_std_pct = where( ismissing(budget_terms_std_pct), 0., budget_terms_std_pct )  ; change missing values to 0 for printing
  delete(denom)

  ;--------------------------------------
  ; Relative diff of multi-year mean
  ;--------------------------------------
  denom = where( budget_terms_mean(:,:,:,0).ne.0, budget_terms_mean(:,:,:,0), FillValue)
  denom@_FillValue = FillValue
  reldiff_of_mean = budget_terms_mean(:,:,:,1)/denom - 1.
  reldiff_of_mean = where( ismissing(reldiff_of_mean), 0., reldiff_of_mean )  ; change missing values to 0 for printing
  delete(denom)

  ;******************************************************
  ; Print table(s) of budget by species
  ;******************************************************
  species_table = new( nrow+2, "string")  ; the 2 extra rows are title rows
  species_table(:) = ""

  if (.not.l_debug) then
     nspectmp = 0                     ; only print table for the first species
  else
     nspectmp = dimsizes(species)-1   ; print tables for all species
  end if

  do isp = 0,nspectmp

     print("")
     print("=========================================")
     print(" Budget table for "+species_longname(isp))
     print("=========================================")
     print("")

     print("\middlehline")

     ; first column of table

     species_table(0)  = "\multirow{ 2}{*}{\bf "+species_longname(isp)+" budget} &"
     species_table(1)  = " & "
     species_table(2:) = row_titles + " & "

     ; header row 0: different regions

     irwh = 0
     do iregion = 0,nregion-1
        species_table(irwh) = species_table(irwh)+ " \multicolumn{3}{c}{{\bf "+regionlbl(iregion)+"}} &&"
     end do

     ; header row 1: 3 columns of numbers in each region

     irwh = irwh+1
     do iregion = 0,nregion-1
        species_table(irwh) = species_table(irwh)+ " "+explbl4table(0)+" & "+explbl4table(1)+" & Rel. diff. &&"
     end do

     ; numbers 

     irwc0 = irwh+1    ; first row that contains numbers

     do ii = 0,nrow-1
     do iregion = 0,nregion-1

        if (ii.lt.nrow-2) then
           floatFmt = "%-.0f"     ; Process rates are large. Show them as integers. 
        else
           floatFmt = "%-.1f"     ; Burden and lifetime are small numbers. Show one digit after decimal point.
        end if

        pcntgFmt = "%-.0f"        ; Show percentages as integers

        str0a = sprintf_PlsMns_LaTex( floatFmt, budget_terms_mean(isp,iregion,ii,0) )
        str1a = sprintf_PlsMns_LaTex( floatFmt, budget_terms_mean(isp,iregion,ii,1) )
        str2  = sprintf_PlsMns_LaTex( pcntgFmt, reldiff_of_mean(isp,iregion,ii)*100 )

        str0b = sprintf( pcntgFmt, budget_terms_std_pct(isp,iregion,ii,0) )
        str1b = sprintf( pcntgFmt, budget_terms_std_pct(isp,iregion,ii,1) )

        species_table(ii+irwc0) = species_table(ii+irwc0) + str0a+" $(\pm"+str0b+"\%)$" \
                                                    +" & "+ str1a+" $(\pm"+str1b+"\%)$" \
                                                    +" & " +str2+"\% &&"
     end do
     end do
     species_table(:) = species_table(:) + "\\"

     print(""+species_table)
     print("")
  end do

