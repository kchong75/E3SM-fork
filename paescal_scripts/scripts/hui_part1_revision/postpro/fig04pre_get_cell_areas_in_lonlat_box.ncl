  load "./utils/sub_lonlat_utils.ncl"

  casename = (/"cflx_1_L72_F20TRC5-CMIP6_6h", \
               "cflx_2_L72_F20TRC5-CMIP6_6h"  /)

  histFilePath = "/pscratch/sd/h/huiwan/cflx/v1_papers/v1_runs/"+casename+"/run/"
  year = 2000

  lverbose = False

  ;=======================================================================================================
  ; Open history file with global output; get cell area and lat/lon info 
  ;=======================================================================================================
  FileName = casename(0)+".cam.h0."+year+"-01.nc"
  File = addfile( histFilePath(0)+FileName,"r")

  area_global = File->area
   lon_global = File->lon
   lat_global = File->lat

  print("")
  print("Global output:")
  print("ncol = "+dimsizes(lon_global))
  print("lon min,max = "+min(lon_global)+", "+max(lon_global))
  print("lat min,max = "+min(lat_global)+", "+max(lat_global))

  ;=======================================================================================================
  ; Get lonlat domain information from history file
  ;=======================================================================================================
  FileName = casename(0)+".cam.h1.2000-01-29-00000.nc"
  File = addfile( histFilePath(0)+FileName,"r")

  ib = 0
  lonlat_box = read_lonlat_box_info( File, ib )

  varname = "lat"+lonlat_box@varname_suffix
  lat_in_box = File->$varname$
 
  varname = "lon"+lonlat_box@varname_suffix
  lon_in_box = File->$varname$ 

  print("")
  print("History file lonlat box = "+lonlat_box@varname_suffix)
  print("ncol = "+dimsizes(lon_in_box))
  print("lon min,max = "+min(lon_in_box)+", "+max(lon_in_box))
  print("lat min,max = "+min(lat_in_box)+", "+max(lat_in_box))
  print("")

  area_in_box = lon_in_box
  area_in_box@long_name = area_global@long_name
  delete(area_in_box@units)

  tol = 0.01
  ncell_in_box = dimsizes(lon_in_box)
  print("Getting information for "+ncell_in_box+" grid cells...")

  do icell = 0,ncell_in_box-1

     idx_global = ind(     abs( lon_in_box(icell) - lon_global ).lt.tol \
                      .and.abs( lat_in_box(icell) - lat_global ).lt.tol )

     if (.not.all(ismissing(idx_global))) then

        area_in_box(icell) = (/area_global(idx_global)/)

        if (lverbose) then
           print("cell "+icell+", lon = "+lon_in_box(icell)+" vs "+lon_global(idx_global)+ \
                               ", lat = "+lat_in_box(icell)+" vs "+lat_global(idx_global)+ \
                              ", area = "+area_in_box(icell) )
        end if

     else
        print("cell "+icell+" not found in global output! Abort.")
        exit
     end if

  end do

  print("area global: min, max = "+min(area_global)+", "+max(area_global))
  print("area in box: min, max = "+min(area_in_box)+", "+max(area_in_box))

  ;======================================================================================================= 
  ; Write NetCDF file
  ;======================================================================================================= 
  outputFileName = "lonlat"+lonlat_box@varname_suffix+".nc"
  print("Writing output file "+outputFileName)

  system("rm -f "+outputFileName)
  fileOut = addfile(outputFileName,"c")

  varname = "lon"+lonlat_box@varname_suffix
  fileOut->$varname$ = lon_in_box
  print(" lon done")

  varname = "lat"+lonlat_box@varname_suffix
  fileOut->$varname$ = lat_in_box
  print(" lat done")

  varname = "area"+lonlat_box@varname_suffix
  fileOut->$varname$ = area_in_box
  print(" area done")
