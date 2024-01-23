function read_emis_array_and_select_time_and_col( emis_min, emis_max, lonlat_subbox, \
                                                  chkptname, nqoi, QoI_name, varname_suffix, \
                                                  File, itimeS, itimeE, dtime)

local iqoi, varname, emis_in_subbox_2d, emis_in_subbox_1d
begin

  col_idx = lonlat_subbox@col_idx_in_hist_lonlat_box

  ;----------------------------------------------------------------
  ; Read in emission at all grid cells in the LONLAT subbox 
  ;----------------------------------------------------------------
  ; Read in the vertically integrated mixing ratio increment corresponding to emission.
  ; Unit: kg/kg * kg/m^3 * m = kg/m^2

  do iqoi = 0,nqoi-1
     varname = QoI_name(iqoi)+"_v_"+ chkptname+"_inc"+varname_suffix
     if (iqoi.eq.0) then ; read the first QoI
        emis_in_subbox_2d = File[:]->$varname$(itimeS:itimeE,col_idx)
     else ; accumulate values 
        emis_in_subbox_2d = emis_in_subbox_2d + File[:]->$varname$(itimeS:itimeE,col_idx)
     end if
  end do

  ; merge time and col dimensions and convert increment to tendency 

  emis_in_subbox_1d = ndtooned(emis_in_subbox_2d)/dtime

  ;----------------------------------------------------------------------------
  ; Select time steps and grid columns in the (emis_min, emis_max] range
  ;----------------------------------------------------------------------------

  coltime_idx = ind( (emis_in_subbox_1d.gt.emis_min) .and. \
                     (emis_in_subbox_1d.le.emis_max)       )

  print("")
  print(" Selecting from  "+(dimsizes(col_idx))+" cols " )
  if (.not.all(ismissing(coltime_idx))) then
     print(" # of grid cols and timesteps selected = "+dimsizes(coltime_idx))
  else
     print(" # of grid cols and timesteps selected = 0! Abort.")
     exit
  end if

  emis_selected_1d = emis_in_subbox_1d(coltime_idx)
  emis_selected_1d@npts = dimsizes(coltime_idx)
  emis_selected_1d@coltime_idx_in_lonlat_subbox = coltime_idx

  var_in_subbox_1d = ndtooned(conform( emis_in_subbox_2d, lonlat_subbox@area, 1 ))
  emis_selected_1d@area = var_in_subbox_1d( coltime_idx )

  var_in_subbox_1d = ndtooned(conform( emis_in_subbox_2d, lonlat_subbox@lon, 1 ))
  emis_selected_1d@lon = var_in_subbox_1d( coltime_idx )

  var_in_subbox_1d = ndtooned(conform( emis_in_subbox_2d, lonlat_subbox@lat, 1 ))
  emis_selected_1d@lat = var_in_subbox_1d( coltime_idx )

  delete(emis_in_subbox_2d)
  delete(emis_in_subbox_1d)
  delete( var_in_subbox_1d)

  print(" Emission  min max = "+(min(emis_selected_1d))+", "+max(emis_selected_1d))
  print(" Cell area min max = "+(min(emis_selected_1d@area))+", "+max(emis_selected_1d@area))

  return(emis_selected_1d)

end
