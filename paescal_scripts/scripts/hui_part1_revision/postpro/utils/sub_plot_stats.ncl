;----
undef("make_boxplot")
procedure make_boxplot( stats_in, scale_fac, chkpt_label, ilevS, ilevE, \
                        plotFileName, polycolor, \
                        trXMinF, trXMaxF, x_axis_scale, tiXAxisString, PanelMainString )

local l_debug, nlev, nchkpt, panel, box1dum, box2dum, linedum, meandum, zerodum, \
      wks, res, FontHeightF, ichkpt, ilev
begin

  l_debug = True
  l_debug = False

  ;-----------------------------
  ; dimension sizes
  ;-----------------------------
  nlev = ilevE - ilevS + 1
  nchkpt = dimsizes(chkpt_label)

  ;-----------------------------------------
  ; Arrays and tmp variables for plotting
  ;-----------------------------------------
  panel   = new(   nchkpt,         "graphic" )
  box1dum = new( (/nchkpt,nlev/),  "graphic" )   ; filled boxes
  box2dum = new( (/nchkpt,nlev,4/),"graphic" )   ; unfilled boxes
  linedum = new( (/nchkpt,nlev,3/),"graphic" )   ; lines and whiskers
  meandum = new( (/nchkpt,nlev/),  "graphic" )   ; mean value
  zerodum = new( (/nchkpt/),       "graphic" )   ; zero lines

  plotFmt = "pdf"
  wks  = gsn_open_wks(plotFmt,plotFileName)
  system("chmod a+r "+plotFileName+"."+plotFmt)
  gsn_define_colormap(wks,"StepSeq25")

  res = True

  FontHeightF = 0.018
  res@tiMainFontHeightF    = FontHeightF*1.
  res@tiXAxisFontHeightF   = FontHeightF
  res@tiYAxisFontHeightF   = FontHeightF
  res@tmXBLabelFontHeightF = FontHeightF*0.8
  res@tmYLLabelFontHeightF = FontHeightF*0.9
  res@gsnStringFontHeightF = FontHeightF
 ;res@tiMainFont = 22 ;helvetica-bold

  res@tmYLTickEndF = (ilevE+1)
  res@trYReverse = True

  res@gsnDraw    = False   ; Don't draw
  res@gsnFrame   = False   ; Don't advance frame
  res@vpWidthF   = 0.4     ; 
  res@vpHeightF  = 0.3
  
  res@trXMinF = trXMinF /x_axis_scale
  res@trXMaxF = trXMaxF /x_axis_scale
  if (abs(res@trXMinF).lt.abs(res@trXMaxF)*0.1) then
     res@trXMinF = -0.1*abs(res@trXMaxF)
  end if
 ;res@tmXBFormat     = "@!;*^se"
  res@tmXBFormat     = "0f"

  res@trYMinF        = ilevS + 0.1 ;+ 5
  res@trYMaxF        = ilevE + 1.8 
  res@tmYLMode       = "Manual"
  res@tmYLTickStartF = ilevS +1
  res@tmYLTickEndF   = ilevE+1
  res@tmYLFormat     = "f"

  res@tmYLTickSpacingF  = 1
  res@gsnTickMarksPointOutward = True
   
 ;res@tmYROn  = False
  res@tmXTOn  = False   ; Turn off top tickmarks

  ;--------------------------------------------------------------------------------------------
  ; Scale the incoming values, so that we can use the same x-axis range for different regions
  ;--------------------------------------------------------------------------------------------
  stats = stats_in * scale_fac /x_axis_scale

  res@tiXAxisString = tiXAxisString
  ;--------------------------------------------------------------------------------------------
  ; Plot statistics in different levels and for different checkpoints
  ;--------------------------------------------------------------------------------------------
  do ichkpt = 0,nchkpt-1

     if (l_debug) then
        print("")
        print(" "+chkpt_label(ichkpt))
        print("")
     end if

     if (ichkpt.eq.0) then 
        res@tiYAxisString = "Model layer index"
     else
        res@tiYAxisString = ""
     end if

     ;-------------------------------------------------------
     ; base plot
     ;-------------------------------------------------------
     res@tiMainString  = chkpt_label(ichkpt)+", x"+scale_fac
    ;res@gsnLeftString = chkpt_label(ichkpt)
    ;res@gsnRightString = " x"+scale_fac

    ;panel(ichkpt) = gsn_blank_plot(wks,res)
     xx00 =  (/res@trXMinF,res@trXMinF/)
     yy00 =  (/res@trYMinF,res@trYMinF/)
     panel(ichkpt) = gsn_csm_xy(wks,xx00,yy00,res)
     zerodum(ichkpt) = gsn_add_polyline( wks,panel(ichkpt),(/0.,0./),(/res@trYMinF,res@trYMaxF/),False)

     ;-----------------
     ; Draw boxes 
     ;-----------------
     boxLineThicknessF = 1.5

     plyres = True
     plyres@gsLineThicknessF = boxLineThicknessF 
   
     mkres = True
     mkres@gsMarkerIndex = 16

     hgttmp = (res@trYMaxF - res@trYMinF)/nlev/1.4
     boxhalfhgt = min( (/0.35,hgttmp/) )

     yoffset = 0.0001 ; otherwise will have problem drawing horizontal lines

     do ilev = ilevS,ilevE 

        if (l_debug) then
           print("")
           print("lev = "+ilev)
        end if

        ; y-axis is model layer index

        yy0 = ilev-boxhalfhgt
        yy1 = ilev+boxhalfhgt
        yy = 1 + (/ yy0, yy0+yoffset, yy1, yy1+yoffset, yy0/)

        ;--------------------------------------------------------
        ; get statistics
        ;--------------------------------------------------------
        kk = ilev - ilevS

        zavg = stats(ichkpt,kk,0) 
        zmed = stats(ichkpt,kk,8)
        zsig = stats(ichkpt,kk,1) 
        zmin = stats(ichkpt,kk,2) 
        zmax = stats(ichkpt,kk,14)

        zqlo = stats(ichkpt,kk,6)  ; lower  quartile
        zqhi = stats(ichkpt,kk,10) ; higher quartile

        zolo = stats(ichkpt,kk,4)  ; lower  octile
        zohi = stats(ichkpt,kk,12) ; higher octile

        zslo = stats(ichkpt,kk,5)  ; lower  sextile
        zshi = stats(ichkpt,kk,11) ; higher sextile

        zdlo = stats(ichkpt,kk,3)  ; lower  decitile
        zdhi = stats(ichkpt,kk,13) ; higher decitile

        if (l_debug) then
           print("lev = "+ilev+", min, max = "+zmin/scale_fac+", "+zmax/scale_fac)
           print("lev = "+ilev+", avg = "+zavg/scale_fac+", median = "+zmed/scale_fac)
        end if

        ;------------
        ; Filled box 
        ;------------
        xx = (/ zqlo, zqhi, zqhi, zqlo, zqlo/)

        plyres@gsFillColor = polycolor(ichkpt)+1
        if (zsig.eq.0) then 
           box1dum(ichkpt,ilev-ilevS) = gsn_add_polyline( wks,panel(ichkpt),xx,yy,plyres )
        else
           box1dum(ichkpt,ilev-ilevS) = gsn_add_polygon ( wks,panel(ichkpt),xx,yy,plyres )
        end if

        if (l_debug) then
           print("lev = "+ilev+", filled box: " +zqlo/scale_fac+", "+zqhi/scale_fac)
        end if

        ;--------------
        ; Unfilled box
        ;--------------
        xx = (/ zslo, zshi, zshi, zslo, zslo/)

        plyres@gsLineThicknessF = boxLineThicknessF
        plyres@gsLineColor =  polycolor(ichkpt)-1

        if (zsig.eq.0) then 
           box2dum(ichkpt,ilev-ilevS,0) = gsn_add_polyline( wks,panel(ichkpt),xx,yy,plyres )
        else
          do ii = 0,3
           box2dum(ichkpt,ilev-ilevS,ii) = gsn_add_polyline ( wks,panel(ichkpt),xx(ii:ii+1),yy(ii:ii+1),plyres )
          end do
        end if

        if (l_debug) then
           print("lev = "+ilev+", unfilled box:" +zslo/scale_fac+", "+zshi/scale_fac)
        end if

        ;-------------------
        ; Line and whiskers 
        ;-------------------
        plyres@gsLineThicknessF = boxLineThicknessF
        plyres@gsLineColor =  polycolor(ichkpt)-1

       ;linedum(ichkpt,ilev-ilevS,0) = gsn_add_polyline( wks,panel(ichkpt),(/zmin,zmax/),(/ilev+1,ilev+1+yoffset/),plyres)

       ;zdlo = stats(ichkpt,kk,24)
       ;zdhi = stats(ichkpt,kk,25)
        linedum(ichkpt,ilev-ilevS,0) = gsn_add_polyline( wks,panel(ichkpt),(/zdlo,zdhi/),(/ilev+1,ilev+1+yoffset/),plyres)

        linedum(ichkpt,ilev-ilevS,1) = gsn_add_polyline( wks,panel(ichkpt),(/zdlo,zdlo/),(/ilev+1-0.7*boxhalfhgt,ilev+1+0.7*boxhalfhgt/),plyres)
        linedum(ichkpt,ilev-ilevS,2) = gsn_add_polyline( wks,panel(ichkpt),(/zdhi,zdhi/),(/ilev+1-0.7*boxhalfhgt,ilev+1+0.7*boxhalfhgt/),plyres)

        if (l_debug) then
           print("lev = "+ilev+", whiskers:" +zdlo/scale_fac+", "+zdhi/scale_fac)
          ;print("lev = "+ilev+", whiskers:" +zdlo+", "+zdhi)
        end if

        ;------
        ; Mean
        ;------
        mkres@gsMarkerColor = 1
        mkres@gsMarkerThicknessF = 1.
        mkres@gsMarkerIndex = 16
        mkres@gsMarkerSizeF = 0.002

       ;meandum(ichkpt,ilev-ilevS) = gsn_add_polymarker( wks,panel(ichkpt),zavg,ilev+1,mkres)
        meandum(ichkpt,ilev-ilevS) = gsn_add_polymarker( wks,panel(ichkpt),zmed,ilev+1,mkres)

     end do 

  end do ; chkpt loop (ichkpt)


  print("")
  print("======================")
  print("")

  resP = True
  resP@gsnMaximize = True

 ;resP@gsnPanelYWhiteSpacePercent = 10 
  resP@gsnPanelXWhiteSpacePercent = 1

  resP@gsnPanelMainFontHeightF = res@tiMainFontHeightF *0.5
  resP@gsnPanelMainString      = PanelMainString
  resP@gsnPanelMainFont        = 22 ;helvetica-bold

  gsn_panel(wks,panel,(/1,nchkpt/),resP)

  delete(panel)
  delete(wks)

end
