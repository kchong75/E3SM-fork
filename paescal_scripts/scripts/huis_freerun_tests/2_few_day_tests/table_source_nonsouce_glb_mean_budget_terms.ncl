  ;-------------------------------------
  ; constants used for unit conversion
  ;-------------------------------------
  pi             = 3.14159265358979323846  ; from E3SM
  earth_radius   = 6.37122e6               ; from E3SM
  earth_sfc_area = 4.*pi*earth_radius^2

  kg2Tg  = 1e-9
  sec2yr = 86400.*365
  ;-------------------------------------

  itime = 0
  year  = 2010

  ;----------------------------------------------------
  l_final_paper = True  ; will print less information 

  explbl4table = (/"opt1","opt2","opt3","opt4","opt40","opt41","opt42","opt43","opt44"/) 
  expname = explbl4table+"_F20TRC5-CMIP6" 
  nexp = dimsizes(expname)

  ;------------------------
  FilePath = "/compyfs/wanh895/scidac_scratch/2023_cflx/2023code/few_day_tests/" + expname + "/run/"
 ;FileName = expname+".cam.h0."+year+".ANN.nc"
  FileName = "case.cam.h0.2009-10-03-00000.nc"

  print("")
  print("=====================")
  print(" simulations being compared:")
  print(" "+expname)
  print("=====================")

  ;------
  species_longname = (/"Dust","Sea salt","MOA","BC","POA"/)
  species          = (/ "dst",     "ncl","mom","bc","pom"/)
  nspecies         = dimsizes(species)
  nspecies         = 1 

  ;-----------------------------------
  ; Aerosol processes or states
  ;-----------------------------------
  label      = (/"Sfc emis.","Dry removal", "Turb. dep.","Grav. Setl.","Act./resus.","Wet removal", "Burden"/)
  inc_fac    = (/ 1.,                  -1.,          -1.,          -1.,           1.,           1., 1./sec2yr/)
  QoI_prefix = (/"SF",                  "",           "",           "",           "",           "", "cnd01_"/)
  QoI_suffix = (/"",                 "DDF",        "TBF",        "GVF","_mixnuc1",      "SFWET", "_v_AERWETRM"/)

  nprocess = dimsizes(QoI_suffix)

  ;-----------------------------------
  ; Types of region to look at
  ;-----------------------------------
  regionlbl = (/"Source-region mean","Non-source-region mean","Global mean"/)
  emis_min = (/ 0.,                   -1.,                     -1./)
  emis_max = (/10.,                    0.,                     10./)
  nemisrange     = dimsizes(emis_max)

  ;------
  spatial_avg = new( (/nspecies,nemisrange,nprocess,nexp/), "float" )

  str_table = new( (/nprocess+1,nexp+1/), "string" )
  str_table(0,1:) = explbl4table
  str_table(1:,0) = label

  ;-----------------------------
  ; open data files
  ;-----------------------------
  File = addfiles(FilePath+FileName,"r")

  area = tofloat(File[0]->area)
  ncol = dimsizes(area)

  print("")
  do isp = 0,nspecies-1

     print("")
     print(" SPECIES = "+species(isp))
     print("")

     ; Get the list of tracer names for this species. 
     ; List is saved in variable QoI_name.

     load "./mam4_species_and_modes.ncl"
     nqoi = dimsizes(QoI_name)  

     spatial_avg_of_mode = new( (/nqoi,nprocess,nexp/), "float" )
     spatial_avg_of_mode = 0.

     ;===============================================================================
     ; Loop over all aerosol processes FOR THIS SPECIES
     ;===============================================================================
     do iprocess = 0,nprocess-1
 
         ;-----------------------------------------------
         ; For this process: get the sum over all modes 
         ;-----------------------------------------------
         ; initialize sum of all modes

         var2d = new( (/ncol,nexp/),"float")
         var2d = 0. 

         ; read in data by mode, sum up, and calculate global averages

         do iqoi = 0,nqoi-1
  
           varname = QoI_prefix(iprocess)+QoI_name(iqoi)+QoI_suffix(iprocess)
  
           if (isfilevar(File[0],varname)) then

              do iexp = 0,nexp-1
                 var2d_of_mode = File[iexp]->$varname$(itime,:)*inc_fac(iprocess)
                 spatial_avg_of_mode(iqoi,iprocess,iexp) = sum(var2d_of_mode*area)/sum(area) * earth_sfc_area * kg2Tg * sec2yr 
                 var2d(:,iexp) = var2d(:,iexp) + var2d_of_mode
              end do
           else
              print(varname+" not found. Skip reading.")
           end if
  
           print("Done reading "+varname)
         end do

         ;---------------------
         ; save emissions
         ;---------------------
         if (QoI_prefix(iprocess).eq."SF") then
            emissions_global = var2d
            print("* Emissions saved")
         end if

         ;--------------------------------------------------------------
         ; use emissions to sort grid columns, then calculate averages 
         ;--------------------------------------------------------------
         do iexp = 0,nexp-1
         do ibin = 0,nemisrange-1

            idx = ind( emissions_global(:,iexp).gt.emis_min(ibin) .and. \
                       emissions_global(:,iexp).le.emis_max(ibin) )

            if (.not.(all(ismissing(idx)))) then
               ;print(" found "+dimsizes(idx)+" grid columns in emission range "+emis_min(ibin)+" to "+emis_max(ibin))
                spatial_avg(isp,ibin,iprocess,iexp) = sum(var2d(idx,iexp)*area(idx))/sum(area(idx)) * earth_sfc_area * kg2Tg * sec2yr
            end if
            delete(idx)
         end do
         end do
 
         delete(var2d) 
     end do
     ;===============================================================================
     ; Done reading and calculating process rates for all processes FOR THIS SPECIES
     ;===============================================================================

    ;do ibin = 0,nemisrange-1  ; source, non-source, global
     do ibin = 0,1  ; source region only

        print("")
        print("======")
        print(regionlbl(ibin))
        print("======")

        do iprocess = 0,nprocess-1

           print("")
           print(" process = "+label(iprocess))
           print(explbl4table+": "+spatial_avg(isp,ibin,iprocess,:))
        end do
     end do

exit

   ; print("")
   ; do iprocess = 0,nprocess-1
   ; print("Global mean: "+sprintf("%10.2f",spatial_avg(isp,2,iprocess,0))+" Tg/yr vs "\
   ;                      +sprintf("%10.2f",spatial_avg(isp,2,iprocess,1))+" Tg/yr from "+label(iprocess) )
   ; end do 

   ; print("")
   ; do iprocess = 0,nprocess-1
   ; print("Src region mean: "+sprintf("%10.2f",spatial_avg(isp,0,iprocess,0))+" Tg/yr vs "\
   ;                          +sprintf("%10.2f",spatial_avg(isp,0,iprocess,1))+" Tg/yr from "+label(iprocess) )
   ; end do 

   ; print("")
   ; do iprocess = 0,nprocess-1
   ; print("Non-src region mean: "+sprintf("%10.2f",spatial_avg(isp,1,iprocess,0))+" Tg/yr vs "\
   ;                              +sprintf("%10.2f",spatial_avg(isp,1,iprocess,1))+" Tg/yr from "+label(iprocess) )
   ; end do 


  delete(QoI_name) 
  end do ;species
  print("")
  print("==================================")

