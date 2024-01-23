load "./utils/sub_lonlat_utils.ncl"
load "./utils/sub_sort_emissions.ncl"
load "./utils/sub_plot_emission_pdf.ncl"
load "./utils/sub_cal_stats.ncl"
load "./utils/sub_plot_stats.ncl"
load "./utils/sub_scatter_plots.ncl"

;-----------------------------------------------------------
dtime = 1800. 
cflx_opt = 1  ; the original coupling scheme

nlevEAMList = (/30,72/)          ; # of layers in the vertical grid
ngroup = dimsizes(nlevEAMList)

stateOrTendList = (/"tend","state"/)
nvartype = dimsizes(stateOrTendList)

;----------------------
do ivartype = 0,nvartype-1
 state_or_tend = stateOrTendList(ivartype)

 do igroup = 0,ngroup-1

  nlevEAM = nlevEAMList(igroup) 

  single_casename     = "cflx_"+cflx_opt+"_L"+nlevEAM+"_F20TRC5-CMIP6_6h"
  single_case_chkpt   = (/"CFLX"+cflx_opt,"AERDRYRM","PBCINI", "STCLD","AERWETRM"/)

  single_plotFileSuffix = "_"+single_casename 
  single_tmp_data_dir   = "./data"+"_"+single_casename+"/"

  instOutputFileSuffix = (/".cam.h1.2000-01-29-00000.nc", \
                           ".cam.h1.2000-02-28-00000.nc", \
                           ".cam.h1.2000-03-30-00000.nc"  /)

  single_histFilePath = "/pscratch/sd/h/huiwan/cflx/v1_papers/v1_runs/"+single_casename+"/run/" 
  plotPath = "/global/cfs/projectdirs/m4359/www/huiwan/cflx/2023_rerun_for_gmd/"
  system("mkdir -p "+single_tmp_data_dir)

  ;-----------------------------------------------------------
  l_emission_pdf    = True
  l_emission_pdf    = False

  l_calculate_stats = True
  l_calculate_stats = False

  l_make_box_plots  = False
  l_make_box_plots  = True

  l_make_scatter_plots  = True
  l_make_scatter_plots  = False

  l_debug        = True
  l_debug        = False

  ;--------------------------------------------
  if (nlevEAM.eq.72) then
     ilevS = nlevEAM  - 10 
  else
  if (nlevEAM.eq.30) then
     ilevS = nlevEAM  - 5
  end if
  end if

  ilevE = nlevEAM  -  1

  itimeS = 0
  itimeE = 10  ; use only a small number of time slices in the history files; for testing
  itimeE = -1  ; use all time slices in the history files

  ;-----------------------------------------------------------
  ; Variables to plot are 3D state or tendencies
  ;-----------------------------------------------------------
  if (state_or_tend.eq."tend") then

     inc_suff = "_inc"  ; or "" for state
     inc_fac  = 1./dtime
     field    = "Mixing ratio tendency"
     title_prefix = "Tendency due to"

     plotFileNameMain = "_boxes_tend"
     trXMinF = -2e-9
     trXMaxF = trXMinF*(-1.) 
     panel_scale_fac = (/1., 200., 10000./)

     x_axis_scale = 1e-9
     var_unit_4plot = "10~S~-9~N~ kg kg~S~-1~N~ s~S~-1~N~"

  else
  if (state_or_tend.eq."state") then

     inc_suff = ""  ; or "" for state
     inc_fac  = 1.
     unit     = "kg kg~S~-1~N~"

     field  = "Mixing ratio"
     title_prefix = "Mixing ratio after"

     plotFileNameMain = "_boxes_state"
     trXMinF = -1e-8 
     trXMaxF =  4e-6 
     panel_scale_fac = (/1., 10., 1000./)

     x_axis_scale = 1e-6
     var_unit_4plot = "10~S~-6~N~ kg kg~S~-1~N~"

  else
    print("Value of state_or_tend unrecognized. Abort.")
    exit
  end if
  end if

  emis_unit = "kg m~S~-2~N~ s~S~-1~N~"
  ;-----------------------------------------------------------

  polycolor = (/4,14,19,25,9,4,14,19,25,9/)  ; for different processes

  ;=======================================================================================================
  ; Regions in which the mean vertical profiles will be calculated and plotted
  ;=======================================================================================================
  latmin = (/  5,  5, 10/)
  latmax = (/ 55, 55, 40/)
  lonmin = (/-25,-25,160/)
  lonmax = (/130,130,230/)
  emis_type = (/1,0,0/)     ; 1 = with emission; 0 = no emission

  nregion = dimsizes(emis_type)

  region_label = (/"Source-region","Vicinity-region","Remote-region"/)
  region_panel_lb = (/"","",""/)
  ;-----------------------------------------------------------------------
  ; Emission flux bins for which the budget analysis is done separately
  ;-----------------------------------------------------------------------
 ;emis_bnds   = (/ 0., 2e-8, 2e-6, 5e-5/)  ; The 4 numbers define 3 ranges
  emis_bnds   = (/ 0., 1./)                ; No longer divide; Process all grid columns with emission as one category.
  nemis_bins  = dimsizes(emis_bnds)-1

  ;--------------------------------------------------
  ; Info about the checkpoints
  ;--------------------------------------------------
  label    = (/"Emissions","Dry removal","Resolved transport","Turb-mix+act/resus","Wet removal"/)
  nchkpt   = dimsizes(label)

  ichkpt_emis = ind(label.eq."Emissions")

  print("")
  if (ichkpt_emis.lt.0) then
     print("Error: did not find checkpoint corresponding to emissions - needed for conditional average")
     exit
  else
     print("chkpt corresponding to emissions = "+(ichkpt_emis+1))
  end if

  ;------------------------------------------------------
  ; Aerosol species and the corresponding mixing ratios
  ;------------------------------------------------------
  species  = "dst"
  QoI_name = "cnd01_"+(/"dst_a1","dst_a3"/)
  nqoi     = dimsizes(QoI_name)

  ;=======================================================================================================
  ; Specify a list of high-freq output files to be processed
  ;=======================================================================================================
  FileName = single_casename + instOutputFileSuffix

  File0 = addfile( single_histFilePath+FileName(0),"r")  ; the first file; used for reading coordinate info

if (l_emission_pdf.or.l_calculate_stats.or.l_make_scatter_plots) then

  File  = addfiles(single_histFilePath+FileName,   "r")  ; all files. Will read fields from these.

  ;=======================================================================================================
  ; Get dimension information from history files
  ;=======================================================================================================
  ; lev
  ;------------------
  nlev  = ilevE - ilevS + 1

  ;------------------
  ; time
  ;------------------
  ntime_in_hist = dimsizes(File[:]->time)
  if (itimeE.lt.0) 
     itimeE = ntime_in_hist - 1
  end if
  ntime = itimeE - itimeS + 1
  print("")
  print("Will process "+ntime+" time slices")

  ;------------------
  ; lat/lon
  ;------------------
  ibox_in_hist = 0
  hist_lonlat_box = read_lonlat_box_info( File0, ibox_in_hist )
   varname_suffix = hist_lonlat_box@varname_suffix

  ;------------------
  ; cell area
  ;------------------
  areaFile = addfile("./lonlat"+varname_suffix+".nc","r")
  hist_lonlat_box@area = read_lonlat_box_cellarea( areaFile, hist_lonlat_box )

end if

  ;================================================================================
  ; Process different regions (source, source-vicinity, and remote).
  ; These regions are defined based on their geographical locations (lat/lon)
  ; and instantaneous emission fluxes).
  ;================================================================================
  print("")
  print("nregion = "+nregion)
  do ir = 0,nregion -1

     print("")
     print("~~~~~~~~~~~~~~~~~")
     print("iregion = "+ir)
     print("~~~~~~~~~~~~~~~~~")

     if (l_emission_pdf.or.l_calculate_stats.or.l_make_scatter_plots) then
     ;===================================================================================================================
     ; Among the grid cols written to history (i.e., hist_lonlat_box), select the ones that are located in lon/lat sub-box.
     ;===================================================================================================================
     print("")
     print("history file lonlat box = "+hist_lonlat_box@varname_suffix)
     print("        analysis subbox = "+lonmin(ir)+" to "+lonmax(ir)+" lon; "+latmin(ir)+" to "+latmax(ir)+" lat")

     lonlat_subbox = find_portion_of_hist_lonlat_box( hist_lonlat_box, \
                                                      latmin(ir),latmax(ir),lonmin(ir),lonmax(ir))
     end if

     ;===================================================================================================================
     ; Plot emission PDF for source regions
     ;===================================================================================================================
     if (l_emission_pdf .and. (emis_type(ir).eq.1)) then

        ; Select grid cells with emission higher than cutoff

        emission_cutoff = 0.

        emis_min = emission_cutoff
        emis_max = 1.

        emis_selected = read_emis_array_and_select_time_and_col( \
                                                               emis_min, emis_max, lonlat_subbox, \
                                                               single_case_chkpt(ichkpt_emis), nqoi, QoI_name, \
                                                               varname_suffix, File, itimeS, itimeE, dtime) 

        ; Plot emission PDF and calculate contribution of individual bins to total emission  

        infoString = "# of model timesteps = "+ntime+ \
                     "~C~# of grid columns and timesteps with emission > "+emission_cutoff+\
                     " = "+emis_selected@npts
        XAxisString = "Emission flux ("+emis_unit+")"

        plot_emission_pdf( plotPath, single_plotFileSuffix, emis_selected, emission_cutoff, \
                           infoString, XAxisString, emis_bnds )

        delete(emis_selected)
     end if

     ;===================================================================================================================
     ; Statistical distributions of state/tendency profiles
     ;===================================================================================================================
     if (emis_type(ir).eq.1) then ; source regions 

        emis_bnds_4plot = emis_bnds
        nemis_rng_4plot = dimsizes(emis_bnds) -1

        scale_fac_4plot = panel_scale_fac(0:nemis_rng_4plot-1)

     else
     if (emis_type(ir).eq.0) then ; vicinity or remote regions

        emis_bnds_4plot = (/-1.,0./)
        nemis_rng_4plot = 1

        scale_fac_4plot = (/panel_scale_fac(dimsizes(emis_bnds)-2+ir)/) 

     else
        print(" emis_type = "+emis_type(ir)+" not recognized. Abort.")
        exit
     end if
     end if

     ;------------------------------------------------------------------------------------
     ; Process different regions and different emission ranges in the source regions
     ;------------------------------------------------------------------------------------
     do ie = 0,nemis_rng_4plot-1

         emis_min = emis_bnds_4plot(ie)
         emis_max = emis_bnds_4plot(ie+1)

         statsFileName = "region_"+ir+"_EmisRange_"+ie+"_"+state_or_tend

         print("")
         print(" --- Processing emission range "+ie+": ("+emis_min+", "+emis_max+"] ---")

         ;------------------------------------------------------------
         ; Select data points
         ;------------------------------------------------------------
         if (l_calculate_stats.or.l_make_scatter_plots) then
  
            emis_selected = read_emis_array_and_select_time_and_col( \
                                      emis_min, emis_max, lonlat_subbox, \
                                      single_case_chkpt(ichkpt_emis), nqoi, QoI_name, \
                                      varname_suffix, File, itimeS, itimeE, dtime) 

         end if
         ;------------------------------------------------------------
         ; Calculate statistics
         ;------------------------------------------------------------
         if (l_calculate_stats) then

            ; Calculate statistics for each chkpt and level
 
            stats = calculate_statistics( single_case_chkpt, QoI_name, inc_suff, inc_fac, \
                                          File, varname_suffix, itimeS, itimeE, ilevS, ilevE, \
                                          lonlat_subbox, emis_selected )

            ; Write results to disk

            system("rm "+single_tmp_data_dir+statsFileName+".nc")
            statsFile = addfile(single_tmp_data_dir+statsFileName+".nc","c")
            statsFile->stats = stats

            delete(stats)
         end if

         ;------------------------------------------------------------
         ; Plot statistics
         ;------------------------------------------------------------
         if (l_make_box_plots) then

            ; Read stats from file

            statsFile = addfile(single_tmp_data_dir+statsFileName+".nc","r")
            stats     = statsFile->stats

            ; Set figure file name, title, etc. 

            plotFileName = plotPath+ "boxplot_"+species+"_"+statsFileName+"_levE"+(ilevE+1)+single_plotFileSuffix
            pageTitle    = region_label(ir)+" x"+scale_fac_4plot(ie)

            chkpt_label = title_prefix+" "+str_lower(label)
            tiXAxisString = field+" ("+var_unit_4plot+")"

            PanelMainString = region_panel_lb(ir)+" "+region_label(ir)+" "+str_lower(field)+" statistics"
            if (emis_type(ir).eq.1) then
               PanelMainString = PanelMainString +" (emission flux range: "+emis_min+" to "+emis_max+" "+emis_unit+")"
                             ; + ", # of data pts per layer: "+
            end if

            ; Make boxplots

            make_boxplot( stats, scale_fac_4plot(ie), chkpt_label, ilevS, ilevE, \
                          plotFileName, polycolor, \
                          trXMinF, trXMaxF, x_axis_scale, tiXAxisString, PanelMainString )
            
            delete(stats)
         end if

         ;------------------------------------------------------------
         ;------------------------------------------------------------
         if (l_make_scatter_plots) then

            if (ir.eq.0 .and. ie.eq.1) then
            plotFileName = plotPath+ "scatter_plot_"+species+"_region_"+ir+"_EmisRange_"+ie+"_"+single_plotFileSuffix

            make_scatter_plots(  QoI_name, File, varname_suffix, itimeS, itimeE, ilevE, dtime, cflx_opt, \
                               lonlat_subbox, emis_selected, plotFileName )

            end if

         end if
         ;------------------------------------------------------------
         ; Clean up
         ;------------------------------------------------------------
         if (l_calculate_stats.or.l_make_scatter_plots) then
            delete(emis_selected)
         end if

     end do ; emission ranges

     delete(emis_bnds_4plot)
     delete(scale_fac_4plot)
  end do ; regions
  ;===============

 end do ; experiment group loop (L30 or L72)
end do  ; state or tendency loop
