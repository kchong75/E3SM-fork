procedure plot_emission_pdf( plotPath,plotFileSuffix, emis_selected, emission_cutoff, \
                             infoString, XAxisString, emis_bnds )

local nemis_pts,nemis_bins,wks,pdfres,nbins,emission_pdf,pdfmax,xx,yy, \
      emis_total, emis_frac, cass_count, case_frac, emis_total, im, refdum, txtdum, \
      lnres, txres, txtdum1
begin

  nemis_pts = emis_selected@npts 
  nemis_bins = dimsizes(emis_bnds)-1

  ;----------------------
  ; Calculate 1D PDF
  ;----------------------
  nbins = 100
  emission_pdf = pdfx(log10(emis_selected), nbins, False)  
  pdfmax       = max(emission_pdf)

  ;----------------------
  ; Plot the 1D PDF
  ;----------------------
  plotFile = plotPath+"/emission_PDF"+plotFileSuffix
  plotFmt  = "pdf"
  wks = gsn_open_wks (plotFmt,plotPath+"/emission_PDF"+plotFileSuffix)
  system("chmod a+r "+plotFile+"."+plotFmt)
  gsn_define_colormap(wks,"StepSeq25")

  pdfres  = True
  pdfres@gsnFrame = False
  pdfres@gsnDraw  = False
 ;pdfres@gsnMaximize       = True    ; Maximize plot in frame
  pdfres@xyLineThicknessF  = 8.
  pdfres@xyLineColor       = 4 
  pdfres@vpHeightF = 0.36
  pdfres@vpWidthF  = 0.6

  FontHeightF = 0.015
  pdfres@tiMainFontHeightF    = FontHeightF*1.1
  pdfres@gsnStringFontHeightF = FontHeightF*0.9
  pdfres@tiXAxisFontHeightF   = FontHeightF
  pdfres@tiYAxisFontHeightF   = FontHeightF
  pdfres@tmXBLabelFontHeightF = FontHeightF*0.9
  pdfres@tmYLLabelFontHeightF = FontHeightF*0.9

  pdfres@gsnXYBarChart            = True              ; Create bar plot
 ;pdfres@gsnXYBarChartOutlineOnly = True

  pdfres@gsnYRefLine = 0.
  pdfres@gsnAboveYRefLineColor = pdfres@xyLineColor

  pdfres@tmXTOn  = False
  pdfres@tmYROn  = False
  pdfres@tmXTBorderOn  = False
  pdfres@tmYRBorderOn  = False

  pdfres@trXLog = True
  pdfres@trXMaxF = max(emis_bnds)*2 ;1.1
  pdfres@trXMinF = max((/1.e-18,emission_cutoff/))

  pdfres@trYMaxF = pdfmax*1.6
  pdfres@trYMinF = -0.05*pdfres@trYMaxF 

  pdfres@gsnLeftString = infoString
  pdfres@tiXAxisString = XAxisString 
  pdfres@tiYAxisString = "Freq. of occurrence (%)"

  yy = emission_pdf
  xx = 10.^emission_pdf@bin_center

  pdfplot = gsn_csm_xy (wks, xx, yy, pdfres)

  print("")
  print("  - Histogram done.")

  delete(xx)
  delete(yy)
  delete(emission_pdf)

  ;--------------------------------------------------------------------------------
  ; Calculate contribution of different emission ranges to total emission
  ;--------------------------------------------------------------------------------
  floattype  = typeof(emis_selected@area)
  emis_total = new(nemis_bins+1,floattype)
  emis_frac  = new(nemis_bins  ,floattype)
  case_count = new(nemis_bins  ,floattype)
  case_frac  = new(nemis_bins  ,floattype)

  do im = 0,nemis_bins-1

    idx = ind( ( emis_selected.ge.emis_bnds(im  ) ) .and.  \;
               ( emis_selected.lt.emis_bnds(im+1) )        )

    if (.not.all(ismissing(idx))) then
       emis_total(im) = sum( emis_selected(idx) * emis_selected@area(idx) )
       case_count(im) =    dimsizes(idx)
       case_frac (im) = 1.*dimsizes(idx)/nemis_pts
    end if

    delete(idx) 
  end do

  emis_total(nemis_bins) = sum( emis_selected * emis_selected@area )
  emis_frac = emis_total(0:nemis_bins-1)/emis_total(nemis_bins)
    
  print("  - emission fraction (%) = "+emis_frac*100) 

  ;--------------------------------------------------------------------------------
  ; Add the percentage contributions to the histogram plot.
  ;--------------------------------------------------------------------------------
  refdum = new( (/nemis_bins,2/), "graphic")
  txtdum = new( (/nemis_bins,2/), "graphic")

  do im = 0,nemis_bins-1

     lnres = True
     lnres@gsLineColor = "gray60"
     lnres@gsLineThicknessF  = 4.
     lnres@gsLineDashPattern = 14

     refdum(im,0) = gsn_add_polyline(wks, pdfplot, (/emis_bnds(im+1),emis_bnds(im+1)/), \
                                                   (/pdfres@trYMinF,0.95*pdfres@trYMaxF/),lnres)
     refdum(im,1) = gsn_add_polyline(wks, pdfplot, (/max((/pdfres@trXMinF,emis_bnds(im)*1.7/)),emis_bnds(im+1)/), \
                                                          0.95*(/pdfres@trYMaxF,pdfres@trYMaxF/),lnres)

     txres = True
     txres@txFontHeightF = FontHeightF *0.8
     txres@txAngleF = 90.
     txres@txAngleF =  0.

     xx = emis_bnds(im+1) * 0.2
     yy1 = pdfres@trYMaxF * 0.88
     yy2 = pdfres@trYMaxF * 0.80

    ;percentage_format = "%4.1f"
     percentage_format = "%2.0f"
     if (.not.(ismissing(emis_frac(im)))) then
        txres@txFontColor = pdfres@xyLineColor -1
        txt = sprintf(percentage_format,emis_frac(im)*100.)
        txtdum(im,0) = gsn_add_text(wks, pdfplot,txt+"%",xx,yy1,txres)

        txres@txFontColor = 1
        txt = sprintf(percentage_format,case_frac(im)*100.)
        txtdum(im,1) = gsn_add_text(wks, pdfplot,txt+"%",xx,yy2,txres)
     end if


  end do

  txres@txAngleF = 0.
  xx = emis_bnds(1)*1e-5

  txres@txFontColor = pdfres@xyLineColor -1
  txtdum1 = gsn_add_text(wks, pdfplot,"Contribution to emitted mass: ",xx,yy1,txres)

  txres@txFontColor = 1
  txtdum2 = gsn_add_text(wks, pdfplot,"Percentage of samples:",        xx,yy2,txres)

  delete(xx)
  delete(yy1)

  ;-------
  ; Done 
  ;-------
  gsn_panel(wks,pdfplot,(/1,1/),False)
  print("  - Emission PDF plot done")

  ;-------
  ;Try the "gsn_histogram" function. Does not seem to support binning in log scale
  ;-------
  ;res = True
  ;res@gsnHistogramBinIntervals = 10^(fspan(-19,-1,50))
  ;res@gsnHistogramComputePercentages = True
  ;res@gsnHistogramPercentSign = True
  ;plot = gsn_histogram(wks,emission_1d,res)
  ;-------

  print("")
  print("************")
  print("Plot file is "+plotFile+"."+plotFmt)
  print("************")

end
