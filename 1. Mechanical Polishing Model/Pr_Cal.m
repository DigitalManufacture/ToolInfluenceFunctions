function Pr=Pr_Cal(Pmax,r0,r_inspot) %input mm,rad, pa, mm
Pr=Pmax*sqrt(1-r_inspot.^2/r0^2);
end