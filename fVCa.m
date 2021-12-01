function Out = fVCa(CA,cCae,Rg,Temp,Far)
 Out = 10^(3)*((Rg*Temp)/(2*Far))*log(cCae./CA);
end