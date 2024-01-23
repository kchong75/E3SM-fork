begin

 nlev = (/30,72/)
 nlevmax = max(nlev)
 ngrid = dimsizes(nlev)

 lev = new( (/2,nlevmax/), "float" )
 lev(0,0:nlev(0)-1) = (/ \
     3.64346569404006, 7.59481964632869, 14.3566322512925,  \
    24.6122200042009, 38.2682997733355, 54.5954797416925, 72.0124505460262,  \
    87.8212302923203, 103.317126631737, 121.547240763903, 142.994038760662,  \
    168.225079774857, 197.908086702228, 232.828618958592, 273.910816758871,  \
    322.241902351379, 379.100903868675, 445.992574095726, 524.687174707651,  \
    609.778694808483, 691.389430314302, 763.404481112957, 820.858368650079,  \
    859.53476652503, 887.020248919725, 912.644546944648, 936.198398470879,  \
    957.485479535535, 976.325407391414, 992.556095123291/)

 lev(1,:) = (/ \
    0.123825413055617, 0.182829235506849, 0.269948862120937,  \
    0.398581703622882, 0.588509146564753, 0.8689385700405, 1.28299490825493,  \
    1.89435247940639, 2.79702693529326, 4.12983315502484, 5.96844936711808,  \
    8.37740437894177, 11.4737872305346, 15.3339382222169, 19.9963379781671,  \
    25.4446965099616, 31.5932512912852, 38.366283094427, 45.6711979394268,  \
    53.3095614343925, 61.0151816703695, 68.4763902297125, 75.355335897261,  \
    81.9462751247347, 88.9105431495609, 96.466673444937, 104.664967286755,  \
    113.560001421463, 123.210991427466, 133.682179932106, 145.0432655327,  \
    157.369877840079, 170.744079521582, 185.254902229151, 200.998939871104,  \
    218.080997033001, 236.61478658749, 256.723685905259, 278.541558309835,  \
    302.213640121883, 327.89750608917, 355.76412943838, 385.999023919911,  \
    418.803462879263, 454.395813989178, 492.468574022527, 531.639531357463,  \
    570.624903103333, 608.643774421584, 645.319967961572, 680.49804443824,  \
    713.704638320434, 744.474829323195, 772.362837943354, 796.952749718736,  \
    817.868823477236, 835.095170886611, 849.661249110595, 863.17642564934,  \
    876.370648334105, 889.222736096271, 901.71183329022, 913.817504420407,  \
    925.519748147335, 936.798994822704, 947.636158854377, 958.012766371339,  \
    967.911080400775, 977.314063869485, 986.205333161068, 993.756956304222,  \
    998.496439491762 /)
 
 ilev = new( (/2,nlevmax+1/), "float" )
 ilev(1,:) = (/ \
          0.100000003197544, 0.147650822913689, 0.218007648100008, \
          0.321890076141867, 0.475273331103897, 0.701744962025609, \
          1.03613217805539, 1.52985763845446, 2.25884732035832, 3.3352065502282, \
          4.92445975982147, 7.01243897441468, 9.74236978346886, 13.2052046776003, \
          17.4626717668335, 22.5300041895007, 28.3593888304225, 34.827113752148, \
          41.905452436706, 49.4369434421476, 57.1821794266375, 64.8481839141014, \
          72.1045965453237, 78.6060752491983, 85.2864750002711, 92.5346112988506, \
          100.398735591023, 108.931198982486, 118.188803860441, 128.233178994491, \
          139.131180869721, 150.95535019568, 163.784405484479, 177.703753558685, \
          192.806050899617, 209.191828842591, 226.970165223411, 246.25940795157, \
          267.187963858948, 289.895152760722, 314.532127483044, 341.262884695297, \
          370.265374181464, 401.732673658358, 435.874252100169, 472.917375878188, \
          512.019772166866, 551.259290548059, 589.990515658606, 627.297033184562, \
          663.342902738581, 697.653186137899, 729.756090502968, 759.193568143422, \
          785.532107743287, 808.373391694185, 827.364255260287, 842.826086512935, \
          856.496411708255, 869.856439590426, 882.884857077784, 895.560615114759, \
          907.863051465681, 919.771957375133, 931.267538919537, 942.330450725871, \
          952.941866982882, 963.083665759796, 972.738495041755, 981.889632697215, \
          990.52103362492, 996.992878983524, 1000. /)

 ilev(0,:30) = (/ \
          2.25523952394724, 5.03169186413288, 10.1579474285245, \
          18.5553170740604, 30.6691229343414, 45.8674766123295, 63.3234828710556, \
          80.7014182209969, 94.9410423636436, 111.69321089983, 131.401270627975, \
          154.586806893349, 181.863352656364, 213.952820748091, 251.704417169094, \
          296.117216348648, 348.366588354111, 409.83521938324, 482.149928808212, \
          567.22442060709, 652.332969009876, 730.445891618729, 796.363070607185, \
          845.353666692972, 873.715866357088, 900.324631482363, 924.964462406933, \
          947.432334534824, 967.538624536246, 985.112190246582, 1000. /)

  ;-----------------------------------------
  ilev_max    = max(ilev)
  ilev_min    = min(ilev)
  ;-----------------------------------------

  plotFileName = "L"+nlev(1)+"_L"+nlev(0)+"_layers"
  wks  = gsn_open_wks("pdf",plotFileName)
  gsn_define_colormap(wks,"cmocean_gray")

  edgecolors = (/2,2/)
  fillcolors = (/240,240/)

   left_bnd    = -1.1
  right_bnd    =  1.1
  box_left  = (/-0.8,0.2/)
  box_right = (/-0.2,0.8/)

  ;-----------------------------------------
  res = True

  FontHeightF = 0.014
  res@tiMainFontHeightF    = FontHeightF*1.
  res@tiXAxisFontHeightF   = FontHeightF
  res@tiYAxisFontHeightF   = FontHeightF
  res@tmXBLabelFontHeightF = FontHeightF*0.8
  res@tmXTLabelFontHeightF = FontHeightF*0.8
  res@tmYLLabelFontHeightF = FontHeightF*0.9
  res@tmYRLabelFontHeightF = FontHeightF*0.9
  res@gsnStringFontHeightF = FontHeightF
  res@tiMainFont = 22 ;helvetica-bold

  res@tmYLTickEndF = ilev_max ;+0.1
  res@trYReverse = True

  res@gsnDraw    = False   ; Don't draw
  res@gsnFrame   = False   ; Don't advance frame
  res@vpWidthF   = 0.4     ; 
  res@vpHeightF  = 0.6
  
 ;res@trYMinF        = ilev_min
  res@trYMinF        = 907.            ; only plot layers below this pressure value
  res@trYMaxF        = ilev_max +0.8

  res@tmYLFormat    = "f"
  res@tmYLOn        = True
  res@tmYROn        = True
  res@tmYRLabelsOn  = True

 ;res@tmYROn        = False
 ;res@tmYRBorderOn  = False

  res@tmXTOn        = False
 ;res@tmXTBorderOn  = False

  res@tmXBOn           = True
 ;res@tmXBMajorLengthF = 0.
 ;res@tmXBMinorLengthF = 0.
  res@tmXBMode         = "Explicit"
  res@tmXBValues       = 0.5*(box_left+box_right)
  res@tmXBLabels       = "L"+nlev

 ;res@tmXBOn         = False
 ;res@tmXBBorderOn  = False
  ;-------------------------------------------------------
  ; base plot
  ;-------------------------------------------------------
  res@tiYAxisString = "Norminal Pressure (hPa)"
  res@tiMainString  = "Lowest model layers in EAMv1" 

  xx00 =  (/left_bnd,right_bnd/)
  yy00 =  (/res@trYMinF,res@trYMinF/)
  panel = gsn_csm_xy(wks,xx00,yy00,res)

  boxdum = new((/  ngrid,nlevmax/), "graphic" )
  middum = new((/2,ngrid,nlevmax/), "graphic" )
  txtdum = new((/  ngrid,nlevmax/), "graphic" )
  sfcdum = new(    ngrid,           "graphic" )

  do igrid = 0,ngrid-1

     print(" Plotting grid with "+nlev(igrid)+" interfaces")

     plyres = True
     plyres@gsFillColor      = fillcolors(igrid)
     plyres@gsEdgesOn        = True
     plyres@gsEdgeColor      = edgecolors(igrid)
     plyres@gsEdgeThicknessF = 2.

     txres = True
     txres@txFontHeightF = FontHeightF*0.6

     left  =  box_left(igrid)
     right = box_right(igrid)

     do kk = 0,nlev(igrid)-1

        ixx = (/left,right,right,left,left/)
        iyy = (/ilev(igrid,kk+1),ilev(igrid,kk+1),ilev(igrid,kk),ilev(igrid,kk),ilev(igrid,kk+1)/)

        txt = sprinti("%02i",kk+1)
        xx_txt = 0.5*(left+right) 
        yy_txt = lev(igrid,kk)

        boxdum(igrid,kk) = gsn_add_polygon(  wks,panel,ixx,iyy,plyres )
        txtdum(igrid,kk) = gsn_add_text(     wks,panel,txt,xx_txt,yy_txt,txres)

        plyres@gsLineThicknessF = 1
        plyres@gsLineDashPattern = 2
        middum(0,igrid,kk) = gsn_add_polyline( wks,panel,(/left,xx_txt-0.07/),(/lev(igrid,kk),lev(igrid,kk)/),plyres )
        middum(1,igrid,kk) = gsn_add_polyline( wks,panel,(/xx_txt+0.07,right/),(/lev(igrid,kk),lev(igrid,kk)/),plyres )
     end do

     plyres@gsLineDashPattern = 0
     plyres@gsLineThicknessF  = 10 

     left  = box_left(igrid)-0.05
     right = box_right(igrid)+0.05
     sfcdum(igrid) = gsn_add_polyline( wks,panel,(/left,right/),(/ilev_max,ilev_max/),plyres )

  end do ; grid loop

  resP = True
  resP@gsnMaximize = True
  gsn_panel(wks,panel,(/1,1/),resP)

end 
