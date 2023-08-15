  ;----------------------------------------------------
  explbl4table = (/"opt1","opt2","opt3","opt4","opt41","opt42","opt43","opt44"/)
  expname = explbl4table+"_F20TRC5-CMIP6"
  casename = expname
  nexp = dimsizes(expname)

  ;------------------------
  FilePath = "/compyfs/wanh895/scidac_scratch/2023_cflx/2023code/1mon/" + expname + "/run/"
  FileName = expname+".cam.h0.2009-10.nc"

  print("")
  print("=====================")
  print(" simulations being compared:")
  print(" "+expname)
  print("=====================")


  itime     = 0
  dtime     = 1800.
  year      = 2009
  month     = 10
  remap_res = "180x360"

  species_category = 1 ; gases
  species_category = 2 ; aerosols

  if (species_category.eq.1) then
     species_longname = (/"O3","H2O2","H2SO4","SO2","DMS","SOAG"/)
     species          = species_longname
     plotFileName     = "state_cmpr_all_gases"
  else
  if (species_category.eq.2) then
     species_longname = (/"Dust","Sea salt","MOA","BC","POA","Sulfate","SOA"/)
     species          = (/"dst","ncl","mom","bc","pom","so4","soa"/)
     plotFileName     = "state_cmpr_all_aerosols"
  else
     print("Unknown species_category. Abort.")
  end if
  end if

  ;-----------------------------------
  ; Paths to and names of climo files
  ;-----------------------------------
  ; FilePath is specified in settings_tend_profiles_2.ncl
  ; casename is specified in settings_tend_profiles_2.ncl

  FileName = "../remap/"+casename+".cam.h0."+year+"-"+month+"."+remap_res+".nc"

  print(""+FileName)

  ;-------------
  ; Color maps
  ;-------------
  state_colormap = "WhiteBlueGreenYellowRed"
   diff_colormap = "nrl_sirkes"

  ;-----------------------------
  ; For plotting
  ;-----------------------------
  ncol = nexp
  exp1 = ispan(1,nexp,1)
  exp0 = exp1-1

  nrow = dimsizes(species)
  nrow = 1 
  plot = new( (/nrow,ncol/),"graphic" )
  wks  = gsn_open_wks("pdf",plotFileName)

  res = True
  res@gsnFrame = False
  res@gsnDraw  = False

  res@cnFillOn       = True
  res@cnLinesOn      = False
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
  vpHgtF_crs = vpWthF_crs *0.55
  res@vpWidthF  = vpWthF_crs 
  res@vpHeightF = vpHgtF_crs

  res@gsnLeftStringOrthogonalPosF = 0.035
  res@gsnRightStringOrthogonalPosF = 0.035

  ;-----------------------------
  ; Panel lay out
  ;-----------------------------
  col_idx = (/"(a","(b","(c"/)
  col_lbl = explbl4table
  row_lbl = species_longname

  ;-----------------------------
  ; open data file
  ;-----------------------------
  File = addfiles(FilePath+FileName,"r")

  lat = File[0]->lat
  lon = File[0]->lon
  lev = File[0]->lev
  lev@long_name = "Norminal Pressure (hPa)"

;=======================================
; Loop through different species
;=======================================
do irow = 0,nrow-1

  isp = irow
  ;----------------------------
  ; History variables to read 
  ;----------------------------
  load "./mam4_species_and_modes.ncl" ; variable QoI_name becomes a list of tracer names for this species
  nqoi = dimsizes(QoI_name)

  ;----------------------------
  ; Contour levels
  ;----------------------------
        qLevels_ref = (/1e-10,2e-10,5e-10,1e-9,2e-9,5e-9,1e-8,2e-8,5e-8,1e-7/) 
    qDiffLevels_ref = (/-5e-8,-1e-8,-5e-9,-1e-9,-5e-10,5e-10,1e-9,5e-9,1e-8,5e-8/) 

  if (any( species(irow).eq.(/"dst","ncl"/) )) then  
        qLevels =     qLevels_ref 
    qDiffLevels = qDiffLevels_ref 
  end if
  if (any( species(irow).eq.(/"mom","bc","pom","SOAG","SO2","H2O2"/) )) then  
        qLevels =     qLevels_ref * 1e-2
    qDiffLevels = qDiffLevels_ref * 1e-2
  end if
  if (species(irow).eq."so4" .or. species(irow).eq."soa") then  
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

  ;-----------------------------------------------------------------------
  ; Read 3d mixing ratio fields, calculate zonal mean, sum up all modes
  ;-----------------------------------------------------------------------
  inc_fac = 1.

  zonavg_all_exps = new( (/nexp,dimsizes(lev),dimsizes(lat)/),"float" )
  zonavg_all_exps = 0.

  do iqoi = 0,nqoi-1

    varname = QoI_name(iqoi)

    do iexp=0,nexp-1
       zonavg_all_exps(iexp,:,:) = zonavg_all_exps(iexp,:,:) \
                                 + dim_avg( File[iexp]->$varname$(itime,:,:,:) )*inc_fac
    end do

  end do 

  ; calculate relative differences

  array2plot        = zonavg_all_exps
  array2plot(0,:,:) = zonavg_all_exps(1,:,:)

  do iexp=1,nexp-1

  array2plot(iexp,:,:) = ( zonavg_all_exps(iexp  ,:,:)   \
                          -zonavg_all_exps(iexp-1,:,:) ) \
                         / zonavg_all_exps(1,     :,:)
  end do

  array2plot!1 = "lev"
  array2plot!2 = "lat"
  array2plot&lev = lev
  array2plot&lat = lat

  res@tmYRMode = "Automatic"
  res@trYLog   = False
  res@trYMinF  = 100. 
  res@trYMaxF  = 1000.

 ;res@tiXAxisString = "Latitude"

  ;----------------------
  ; plot new results 
  ;----------------------
   res@gsnLeftString        = "~F22~"+row_lbl(irow)

   icol = 0
   res@cnFillPalette        = state_colormap
   res@cnLevelSelectionMode = "ExplicitLevels"   ; set explicit contour levels
   res@cnLevels             = qLevels 
  ;res@gsnLeftString        = "~F22~"+col_idx(icol)+(irow+1)+") "+row_lbl(irow)
   res@gsnRightString       = col_lbl(1) 
   plot(irow,icol) = gsn_csm_pres_hgt(wks,array2plot(0,:,:),res)    
   delete(res@cnLevels)

   res@cnFillPalette  = diff_colormap 
  ;-------------------------
  ; plot diff
  ;-------------------------
   res@cnLevelSelectionMode = "ExplicitLevels"   ; set explicit contour levels
  ;res@cnLevels             = qDiffLevels 
   res@cnLevels             = (/-1.5,-1,-0.5,-0.3,-0.1,0.1,0.3,0.5,1,1.5/)

   do icol = 1,nexp-1

  ;res@gsnLeftString        = "~F22~"+col_idx(icol)+(irow+1)+") "+row_lbl(irow)
   res@gsnRightString       = col_lbl(icol) +" minus " + col_lbl(icol-1)

   plot(irow,icol) = gsn_csm_pres_hgt(wks,array2plot(icol,:,:),res)    

   end do     
   delete(res@cnLevels)

delete(QoI_name)
end do

  resP = True
  resP@gsnMaximize = True
  resP@gsnPanelYWhiteSpacePercent = 5
  resP@gsnPanelXWhiteSpacePercent = 3 
 ;gsn_panel(wks,ndtooned(plot),(/nrow,ncol/),resP)
  gsn_panel(wks,ndtooned(plot),(/4,2/),resP)

