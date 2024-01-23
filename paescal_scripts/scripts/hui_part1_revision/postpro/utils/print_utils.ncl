function sprintf_PlsMns_LaTex(  fmtIn:string     \
                              , numberIn:numeric )
;-------------------------------------------------------------------------
local sign
begin

  strOut = sprintf(fmtIn,abs(numberIn))
  if (numberIn.gt.0.) then
     sign = "$+$"
  else
     sign = "$-$"
  end if
     
  strOut = sign + strOut 
  return(strOut)
end


