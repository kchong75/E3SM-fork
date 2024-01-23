
function calculate_statistics( chkpt, QoI_name, inc_suff, inc_fac, \
                               File, varname_suffix, itimeS, itimeE, ilevS, ilevE, \
                               lonlat_subbox, emis_selected )

local nchkpt,ichkpt,nqoi,iqoi,nlev,kk,varname,col_idx,var_in_subbox,coltime_idx, \
      var1d_in_subbox,var_kk_selected,xx,xx_stats
begin

    ;--------------------------------------
    ; Dimension sizes
    ;--------------------------------------
    nchkpt = dimsizes(chkpt)
    nqoi   = dimsizes(QoI_name)
    nlev   = ilevE - ilevS + 1

    ;--------------------------------------
    ; Output array
    ;--------------------------------------
    ; Find how many statistics the function stat_dispersion will provide

    xx = random_uniform(0.,1.,100)
    xx_stats = stat_dispersion(xx,False)
    nstats = dimsizes(xx_stats)
    delete([/xx,xx_stats/])

    ; Create output array; name dimensions for NetCDF output

    stats   = new( (/nchkpt,nlev,nstats/),typeof(emis_selected) )
    stats!0 = "chkpt"
    stats!1 = "lev"

    ;------------------------------------------------------------
    ; Process each chkpt separately
    ;------------------------------------------------------------
    col_idx = lonlat_subbox@col_idx_in_hist_lonlat_box

    do ichkpt = 0,nchkpt-1

     print("")
     print(" ++ Reading data for checkpoint "+chkpt(ichkpt))

     ;------------------------------------------------------------------------------------
     ; read in tendency or mixing ratio in the analysis latlon box and sum up all modes
     ;------------------------------------------------------------------------------------
     do iqoi = 0,nqoi-1
        varname = QoI_name(iqoi)+"_"+ chkpt(ichkpt)+inc_suff+varname_suffix
        if (iqoi.eq.0) then ; read the first QoI
           var_in_subbox = File[:]->$varname$(itimeS:itimeE,ilevS:ilevE,col_idx)*inc_fac
           print("    data reading done for "+varname)
        else ; accumulate values 
           var_in_subbox = var_in_subbox + File[:]->$varname$(itimeS:itimeE,ilevS:ilevE,col_idx)*inc_fac
           print("    data reading done for "+varname)
        end if
     end do
     print("    min = "+min(var_in_subbox)+", max = "+max(var_in_subbox)+" in lonlat_subbox")

     ;------------------------------------------------------------------------------------
     ; Select timesteps and cols according to emis_selected 
     ;------------------------------------------------------------------------------------
     coltime_idx = emis_selected@coltime_idx_in_lonlat_subbox

     ;------------------------------------------------------------------------------------
     ; Calculate statistics 
     ;------------------------------------------------------------------------------------
     print(" ++ Obtaining statistics for checkpoint "+chkpt(ichkpt))
     do kk = 0,nlev-1
        ;print("    lev "+(kk+ilevS))
        var1d_in_subbox = ndtooned( var_in_subbox(:,kk,:) )
        var_kk_selected = var1d_in_subbox( coltime_idx )
        stats(ichkpt,kk,:) = stat_dispersion( var_kk_selected, False )
     end do 

    end do ; checkpoints
    ;--------------------

    stats!2 = "statistics"
    stats@url = "https://www.ncl.ucar.edu/Document/Functions/Contributed/stat_dispersion.shtml"

   return(stats)

end


