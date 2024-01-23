undef("cell_is_in_latlon_box")
function cell_is_in_latlon_box( lat:numeric, lon:numeric, \
                                latmin:numeric, latmax:numeric, \
                                lonmin:numeric, lonmax:numeric )
;-------------------------------------------------------------------------
local cell_is_in_lat_range, cell_is_in_lon_range, lon_minus, lon_plus
begin

  cell_is_in_lat_range =  (lat.ge.latmin).and.(lat.le.latmax)

  lon_minus = lon + 360.
  lon_plus  = lon - 360.

  cell_is_in_lon_range = ( (lon      .ge.lonmin).and.(lon      .le.lonmax) ) .or. \
                         ( (lon_minus.ge.lonmin).and.(lon_minus.le.lonmax) ) .or. \ 
                         ( (lon_plus .ge.lonmin).and.(lon_plus .le.lonmax) )

  flag = where( cell_is_in_lat_range.and.cell_is_in_lon_range, 1, 0 )

  return(flag)
end


undef("read_lonlat_box_info")
function read_lonlat_box_info( aFile, ib )
;-------------------------------------------------------------------------
local file_dim_names,lonlat_dim_names,lonlat_box_label,nlonlat_box,lat_varname,lon_varname,varname_suffix
begin

  ; Find the list of lonlat boxes in the history file.

  file_dim_names   = getfiledimnames( aFile )                      ; get all dimension names
  lonlat_dim_names = str_match( file_dim_names, "ncol_" )          ; select the subset starting with "ncol_"

  if (.not.all(ismissing(lonlat_dim_names))) then

     lonlat_box = True

     ; Get the labels; they are needed for determining variable names later in this script

     lonlat_box_label = str_sub_str( lonlat_dim_names, "ncol_", "" )  ; remove "ncol_" from those strings
     varname_suffix = "_"+lonlat_box_label(ib)

     lonlat_box@varname_suffix = varname_suffix

     ; Get lat/lon information for all cells in the lonlat box to be used

     lat_varname = "lat"+varname_suffix
     lon_varname = "lon"+varname_suffix

     lonlat_box@lat = aFile->$lat_varname$
     lonlat_box@lon = aFile->$lon_varname$

  else
     lonlat_box = False
  end if

  return(lonlat_box)
end

undef("read_lonlat_box_cellarea")
function read_lonlat_box_cellarea( areaFile, hist_lonlat_box )
local varname,lontmp,lattmp,tol,cellarea
begin

  ; Verify lat/lon in the history file match the lat/lon in areaFile

  varname = "lon"+hist_lonlat_box@varname_suffix
  lontmp = areaFile->$varname$

  varname = "lat"+hist_lonlat_box@varname_suffix
  lattmp = areaFile->$varname$

  tol = 0.01
  if (any( abs( lontmp - hist_lonlat_box@lon ).gt.tol .or. \
           abs( lattmp - hist_lonlat_box@lat ).gt.tol      )) then
     print("lon/lat info in history file and areaFile do not match! Abort.")
     exit
  end if

  delete(lattmp)
  delete(lontmp)

  ; After verification, read in the cell areas

  varname = "area"+hist_lonlat_box@varname_suffix
  cellarea = areaFile->$varname$
  print("")
  print("Reading of cell area done")

  return(cellarea)

end

;-------------------------------------------------------------------------
undef("find_portion_of_hist_lonlat_box")
function find_portion_of_hist_lonlat_box( hist_lonlat_box,latmin,latmax,lonmin,lonmax)

local col_idx
begin
     ; Look up locations in hist_lonlat_box

     col_idx = ind(1.eq. cell_is_in_latlon_box(hist_lonlat_box@lat,hist_lonlat_box@lon, \
                                               latmin,latmax,lonmin,lonmax )            )

     if (.not.all(ismissing(col_idx))) then
        print(" # of common grid cols = "+dimsizes(col_idx))
        lonlat_portion = True 
        lonlat_portion@col_idx_in_hist_lonlat_box = col_idx
        lonlat_portion@area = hist_lonlat_box@area( col_idx )
        lonlat_portion@lat  = hist_lonlat_box@lat(  col_idx )
        lonlat_portion@lon  = hist_lonlat_box@lon(  col_idx )

     else
        lonlat_portion = False
        print(" # of common grid cols = 0! Abort.")
        exit
     end if

     delete(col_idx)

     return(lonlat_portion)
end
