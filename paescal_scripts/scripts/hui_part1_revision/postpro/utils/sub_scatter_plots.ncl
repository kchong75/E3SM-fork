
procedure make_scatter_plots(  QoI_name, File, varname_suffix, itimeS, itimeE, ilevE, dtime, opt, \
                               lonlat_subbox, emis_selected, plotFileName )

local nchkpt,ichkpt,nqoi,iqoi,nlev,kk,varname,col_idx,var_in_subbox,coltime_idx, \
      var1d_in_subbox,var_kk_selected,xx,xx_stats
begin

    var = (/ (/"PBCINI",   "STCLD_inc"/), \
             (/"CFLX"+opt, "AERDRYRM_inc"/), \
             (/"CFLX"+opt+"_inc", "STCLD_inc"/) /)

    q_scale = 1e-5
    q_unit = "(10~S~-5 ~N~kg kg~S~-1~N~)"

    dqdt_scale = 1e-8
    dqdt_unit = "(10~S~-8~N~ kg kg~S~-1~N~s~S~-1~N~)"

    var_longname = (/ (/ "Mixing ratio "+q_unit, "Tendency due to turb.mix.+act./resus. "+dqdt_unit/), \
                      (/ "Mixing ratio "+q_unit, "Tendency due to dry removal "+dqdt_unit/) , \
                      (/ "Tendency due to emissions "+dqdt_unit,    "Tendency due to turb.mix.+act./resus. "+dqdt_unit/) /)

    title = (/"(a) C vs. q","(b) B vs. q", "(c) C vs. A"/)+" in lowest model layer"
    draw_1to1_line = (/False, False, False/)
    draw_regression_line = (/False, False, False/)

    scale = (/(/1./q_scale,         1/dtime/dqdt_scale/), \
              (/1./q_scale,         1/dtime/dqdt_scale/), \
              (/1/dtime/dqdt_scale, 1/dtime/dqdt_scale/) /)

    color = (/25,13,36/)

    varmin = (/ (/-1e-5/q_scale, -4e-8/dqdt_scale/), \
                (/-1e-5/q_scale, -4e-8/dqdt_scale/), \
                (/-0.3e-8/dqdt_scale, -3.e-8/dqdt_scale/) /)

    varmax = (/ (/16e-5/q_scale,  0.5e-8/dqdt_scale/), \
                (/16e-5/q_scale,  0.5e-8/dqdt_scale/), \
                (/8e-8/dqdt_scale,  0.5e-8/dqdt_scale/) /)
    ;--------------------------------------
    ; Dimension sizes
    ;--------------------------------------
    shape  = dimsizes(var)
    npanel = shape(0)
    nvar_per_panel = shape(1) 
    if (nvar_per_panel.ne.2) then
       print("nvar_per_panel /= 2. Abort.")
       exit
    end if
    delete(shape)

    nqoi   = dimsizes(QoI_name)

    var_selected = new( (/nvar_per_panel+1,emis_selected@npts/),"float" )

    plot = new( npanel, "graphic" )
    dum0 = new( npanel, "graphic" )
    dum1 = new( npanel, "graphic" )
    dumt = new( npanel, "graphic" )

    plotFmt = "pdf"
    wks  = gsn_open_wks(plotFmt,plotFileName)
    system("chmod a+r "+plotFileName+"."+plotFmt)
    gsn_merge_colormaps(wks,"StepSeq25","precip_11lev")

    print("+++ Making scatter plots of process rates.")
    print("+++ Plot file is "+plotFileName+"."+plotFmt)

    ;------------------------------------------------------------
    ; Process each panel separately
    ;------------------------------------------------------------
    col_idx = lonlat_subbox@col_idx_in_hist_lonlat_box

    do ipanel = 0,npanel-1

       do ivar = 0,nvar_per_panel-1

          print("")
          print(" ++ Reading data for "+ var(ipanel,ivar))

          ;------------------------------------------------------------------------------------
          ; read in data in the analysis latlon box and sum up all modes
          ;------------------------------------------------------------------------------------
          do iqoi = 0,nqoi-1
             varname = QoI_name(iqoi)+"_"+ var(ipanel,ivar) +varname_suffix
             if (iqoi.eq.0) then ; read the first QoI
                var_in_subbox = File[:]->$varname$(itimeS:itimeE,ilevE,col_idx)*scale(ipanel,ivar)
                print("    data reading done for "+varname)
             else ; accumulate values 
                var_in_subbox = var_in_subbox + File[:]->$varname$(itimeS:itimeE,ilevE,col_idx)*scale(ipanel,ivar)
                print("    data reading done for "+varname)
             end if
          end do
          ;print("    min = "+min(var_in_subbox)+", max = "+max(var_in_subbox)+" in lonlat_subbox")

          ;------------------------------------------------------------------------------------
          ; Select timesteps and cols according to emis_selected 
          ;------------------------------------------------------------------------------------
          coltime_idx = emis_selected@coltime_idx_in_lonlat_subbox

          var1d_in_subbox = ndtooned( var_in_subbox )
          var_selected(ivar,:) = var1d_in_subbox( coltime_idx )
          print("    min = "+min(var_selected(ivar,:))+", mediam = "+dim_median(var_selected(ivar,:))\
                 +", max = "+max(var_selected(ivar,:))+" in selected samples")

        end do 

        res = True

        res@gsnDraw  = False
        res@gsnFrame = False

        res@xyMarkLineMode = "Markers"
        res@xyMarker       = 1 
        res@xyMarkerSizeF  = 0.02 
        res@tiXAxisString  = var_longname(ipanel,0) ;var(ipanel,0) 
        res@tiYAxisString  = var_longname(ipanel,1) ;var(ipanel,1)
        res@xyMarkerColor  = color(ipanel)

        res@trXMinF = varmin(ipanel,0)
        res@trXMaxF = varmax(ipanel,0)
        res@trYMinF = varmin(ipanel,1)
        res@trYMaxF = varmax(ipanel,1)

        res@gsnXRefLine = 0.
        res@gsnYRefLine = 0.

        res@tmXTOn = False
        res@tmYROn = False

       ;res@tmYLFormat = "@!;*^se"
       ;res@tmXBFormat = "@!;*^se"
        res@tmYLFormat = "f"
        res@tmXBFormat = "f"

        FontHeightF = 0.018*1.2
        res@tiMainFontHeightF    = FontHeightF*1.
        res@tiXAxisFontHeightF   = FontHeightF
        res@tiYAxisFontHeightF   = FontHeightF
        res@tmXBLabelFontHeightF = FontHeightF*0.8
        res@tmYLLabelFontHeightF = FontHeightF*0.9
        res@gsnStringFontHeightF = FontHeightF

        res@gsnCenterString = title(ipanel)

       ;-------------------
       ; line of regression
       ;-------------------
       if (draw_regression_line(ipanel)) then

        res@xyMarkLineModes = (/"Markers","Lines"/)
        res@xyLineThicknesses   = (/1,2/) 
        res@xyLineColor         = "red"
         
        rc = regline_stats( var_selected(0,:), var_selected(1,:) )
        var_selected(2,:) = rc@Yest

        txt = "y = "+sprintf("3.1f",rc@b(0))+" "+sprintf("5.2f",rc@b(1))+"x,  r~S~2 = ~N~"+rc@r2
        print(txt)
        print("r2 = "+rc@r2)

        tres = True
        tres@txFontHeightF = FontHeightF*0.7
        dumt(ipanel) = gsn_add_txt(wks,plot(ipanel), txt,0.8*res@trXMaxF, 0.4*res@trYMinF, txtres )

        end if
       ;-----

        stride = 3 
        plot(ipanel) = gsn_csm_xy(wks,var_selected(0,::stride),var_selected(1:,::stride),res)

      ;-----------------------
      ; Mean values
      ;-----------------------
      ; mres = True
      ; mres@gsMarkerSizeF = res@xyMarkerSizeF * 5
      ; mres@gsMarkerColor = "red"
      ; mres@gsMarkerIndex = 2
      ;
      ; dum0(ipanel) = gsn_add_polymarker(wks,plot(ipanel), dim_avg(var_selected(0,::stride)), \
      ;                                                     dim_avg(var_selected(1,::stride)), mres )
      ;-----------------------

        if (draw_1to1_line(ipanel)) then

            lres = True
            lres@gsLineDashPattern = 5 
            lres@gsLineThicknessF  = 2 
            lres@gsLineColor       = "red" 
            dum1(ipanel) = gsn_add_polyline(wks,plot(ipanel), (/0,-2.*res@trYMinF/), \
                                                              (/0, res@trYMinF/), lres )
        end if

    end do ;panel 
    ;--------------------

    resp = True
    resp@gsnPanelXWhiteSpacePercent  = 8
    gsn_panel(wks,plot,(/1,npanel/),resp)
end

