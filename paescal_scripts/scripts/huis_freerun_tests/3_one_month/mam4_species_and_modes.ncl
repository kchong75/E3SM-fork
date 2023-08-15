  if (species(isp).eq."dst") then  
     QoI_name = (/"dst_a1",         "dst_a3"         /)
     burden_vname = "BURDENDUST"
  else
  if (species(isp).eq."ncl") then  
     QoI_name = (/"ncl_a1","ncl_a2","ncl_a3"         /)
     burden_vname = "BURDENSEASALT"
  else
  if (species(isp).eq."mom") then  
     QoI_name = (/"mom_a1","mom_a2","mom_a3","mom_a4"/)
     burden_vname = "BURDENMOM"
  else
  if (species(isp).eq."soa") then  
     QoI_name = (/"soa_a1","soa_a2","soa_a3"         /)
     burden_vname = "BURDENSOA"
  else
  if (species(isp).eq."so4") then  
     QoI_name = (/"so4_a1","so4_a2","so4_a3"         /)
     burden_vname = "BURDENSO4"
  else
  if (species(isp).eq."bc") then  
     QoI_name = (/ "bc_a1",          "bc_a3", "bc_a4"/)
     burden_vname = "BURDENBC"
  else
  if (species(isp).eq."pom") then  
     QoI_name = (/"pom_a1",         "pom_a3","pom_a4"/)
     burden_vname = "BURDENPOM"

  else
  if (species(isp).eq."num") then  
     QoI_name = (/"num_a1","num_a2","num_a3","num_a4"/)
     burden_vname = ""

  else ; assuming gas

     QoI_name = (/species(irow)/)
     burden_vname = ""

  end if ;1
  end if ;2
  end if ;3
  end if ;4
  end if ;5
  end if ;6
  end if ;7
  end if ;8

  
