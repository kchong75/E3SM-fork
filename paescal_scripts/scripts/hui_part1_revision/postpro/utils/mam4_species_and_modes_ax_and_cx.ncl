  if (species(isp).eq."dst") then  
     QoI_name = (/"dst_a1", "dst_a3", \
                  "dst_c1", "dst_c3" /)
  else
  if (species(isp).eq."ncl") then  
     QoI_name = (/"ncl_a1","ncl_a2","ncl_a3", \
                  "ncl_c1","ncl_c2","ncl_c3" /)
  else
  if (species(isp).eq."mom") then  
     QoI_name = (/"mom_a1","mom_a2","mom_a3","mom_a4", \
                  "mom_c1","mom_c2","mom_c3","mom_c4"/)
  else
  if (species(isp).eq."soa") then  
     QoI_name = (/"soa_a1","soa_a2","soa_a3", \
                  "soa_c1","soa_c2","soa_c3"  /)
  else
  if (species(isp).eq."so4") then  
     QoI_name = (/"so4_a1","so4_a2","so4_a3", \
                  "so4_c1","so4_c2","so4_c3"  /)
  else
  if (species(isp).eq."bc") then  
     QoI_name = (/ "bc_a1",  "bc_a3", "bc_a4", \
                   "bc_c1",  "bc_c3", "bc_c4"  /)
  else
  if (species(isp).eq."pom") then  
     QoI_name = (/"pom_a1", "pom_a3","pom_a4", \
                  "pom_c1", "pom_c3","pom_c4"  /)

  else
  if (species(isp).eq."num") then  
     QoI_name = (/"num_a1","num_a2","num_a3","num_a4", \
                  "num_c1","num_c2","num_c3","num_c4"/)

  else ; assuming gas

     QoI_name = (/species(irow)/)

  end if ;1
  end if ;2
  end if ;3
  end if ;4
  end if ;5
  end if ;6
  end if ;7
  end if ;8

  
