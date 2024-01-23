;----------------------------------------------------------------------
; coneff_18.ncl
;
; Concepts illustrated:
;   - Using functions for cleaner code
;   - Setting contour line thicknesses and patterns
;   - Using "setvalues" to set resource values
;   - Using "getvalues" to retrieve resource values
;   - Drawing partially transparent filled contours
;----------------------------------------------------------------------
;
; These files are loaded by default in NCL V6.2.0 and newer
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;----------------------------------------------------------------------
; Function that retrieves contour levels from a plot and changes the
; the requested contour levels to the requested patterns.
;
; This code is somewhat complicated because we are checking to make
; sure we don't override any line patterns that might have 
; originally been set by the user.
;
; This code is very similar to the set_contour_line_thicknesses code.
;----------------------------------------------------------------------
undef("set_contour_line_patterns")
procedure set_contour_line_patterns(plot,levels_to_change[*]:numeric,\
                                    patterns[*]:integer)
local ii, levels, nlevels, line_patterns, nchange, ncolor, n, \
      mono_pattern, line_pattern, changed_a_level
begin
;---Retrieve the original line patterns (or pattern) used for the plot.
  getvalues plot@contour
    "cnLevels"              : levels
    "cnLineDashPatterns"    : line_patterns
    "cnMonoLineDashPattern" : mono_pattern
    "cnLineDashPattern"     : line_pattern
  end getvalues

;----------------------------------------------------------------------
; If the original contour plot used a single dash pattern for
; all contour lines, then make sure we use that again for 
; any contour lines that are not being changed. Otherwise,
; we assume the user set an array of dash patterns, and
; we'll use those.
;----------------------------------------------------------------------

  if(mono_pattern) then
    line_patterns = line_pattern
  end if

  nlevels  = dimsizes(levels)
  nchange  = dimsizes(levels_to_change)
  npattern = dimsizes(patterns)

  if(npattern.ne.nchange) then
      print("set_contour_line_patterns: error: the contour line patterns must be an")
      print("array of the same length as the number of contour levels to change.")
  end if

;---Array to hold index values of contour levels that need to have a color applied
  changed_a_level = False
  do n=0,nchange-1
    ii := ind(levels.eq.levels_to_change(n))
    if(ismissing(ii(0))) then
      print("set_contour_line_patterns: warning, no contour level equal to " + levels_to_change(n))
    else
      line_patterns(ii(0)) = patterns(n)
      changed_a_level      = True
    end if   
  end do

  if(changed_a_level) then
    setvalues plot@contour
      "cnMonoLineDashPattern" : False             ; allows an array of line patterns to be set
      "cnLineDashPatterns"    : line_patterns
    end setvalues
  end if
end

;----------------------------------------------------------------------
; Function that retrieves contour levels from a plot and changes the
; the requested contour levels to the requested thicknesses.
;
; This code is somewhat complicated because we are checking to make
; sure we don't override any line thicknesses that might have 
; originally been set by the user.
;
; This code is very similar to the set_contour_line_patterns code.
;----------------------------------------------------------------------
undef("set_contour_line_thicknesses")
procedure set_contour_line_thicknesses(plot,levels_to_change[*]:numeric,\
                                       thicknesses[*]:numeric)
local ii, levels, nlevels, line_thicknesses, nchange, nthick, \
      mono_thickness, line_thickness,changed_a_level
begin
;---Retrieve the original line thicknesses (or thickness) used for the plot.
  getvalues plot@contour
    "cnLevels"            : levels
    "cnLineThicknesses"   : line_thicknesses
    "cnMonoLineThickness" : mono_thickness
    "cnLineThicknessF"    : line_thickness
  end getvalues
;----------------------------------------------------------------------
; If the original contour plot used a single dash line thickness
; all contour lines, then make sure we use that again for 
; any contour lines that are not being changed. Otherwise,
; we assume the user set an array of line thicknesses, and
; we'll use those.
;----------------------------------------------------------------------

;---Check user set line thicknesses, and apply new thicknesses if needed
  if(mono_thickness) then
    line_thicknesses = line_thickness
  end if

  nlevels = dimsizes(levels)
  nchange = dimsizes(levels_to_change)
  nthick  = dimsizes(thicknesses)

  if(nthick.ne.nchange) then
      print("set_contour_line_thicknesses: error: the contour line thicknesses must be an")
      print("array of the same length as the number of contour levels to change.")
  end if

;---Array to hold index values of contour levels that need to have a thickness applied
  changed_a_level = False
  do n=0,nchange-1
    ii := ind(levels.eq.levels_to_change(n))
    if(ismissing(ii(0))) then
      print("set_contour_line_thicknesses: warning, no contour level equal to " + levels_to_change(n))
    else
      line_thicknesses(ii(0)) = thicknesses(n)
      changed_a_level     = True
    end if   
  end do

  if(changed_a_level) then
    setvalues plot@contour
      "cnMonoLineThickness" : False             ; allows an array of line thicknesses to be set
      "cnLineThicknesses"   : line_thicknesses
    end setvalues
  end if
end


;;----------------------------------------------------------------------
;; Main code
;;----------------------------------------------------------------------
;begin
;
;;---Read data file and open workstation
;  a = addfile("$NCARG_ROOT/lib/ncarg/data/cdf/uv300.nc","r")
;  u = a->U(1,:,:)  
;;  data = generate_2d_array(15, 20, 0., 80., 0, (/129,129/))
;
;  wks   = gsn_open_wks ("png", "coneff" )  ; send graphics to PNG file
;
;;---Set plot resources
;  res                = True 
;  res@gsnMaximize    = True
;  res@cnFillOn       = True
;  res@mpFillOn       = False
;  res@tiMainString = "Default plot"
;
;;---First draw the default plot
;  plot = gsn_csm_contour_map(wks, u, res)
;
;;---Recreate plot but with contours more opaque
;  res@gsnDraw        = False    ; Don't draw plot or
;  res@gsnFrame       = False    ;   advance frame
;  res@cnFillOpacityF = 0.5
;  res@tiMainString   = "Change line thicknesses and patterns"
;
;  plot = gsn_csm_contour_map(wks, u, res)   
;
;;---Change some line thicknesses and patterns and redraw the plot
;  set_contour_line_thicknesses(plot,(/20,24,28,32/),(/2.,3.5,2.,7.5/))
;  set_contour_line_patterns(plot,(/20,28/),(/2,10/))
;
;  draw(plot)
;  frame(wks)
;
;end

