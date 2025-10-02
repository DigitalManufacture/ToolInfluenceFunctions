dg=0.0015;
Hw=1100*9.8*1e6;
c=3;% relationship between hardness and yield stress
Esp=150e6;
sigw=linspace(0,dg/2,100000);
minA=abs(-dg^3);
flag=1;
for i=1:length(sigw)
    if abs(sigw(i)^3+(9*pi*pi/8*(Hw/c/2)^2/Esp^2-3)*dg*(sigw(i))^2+3*dg^2*sigw(i)-dg^3)<minA
        minA=abs(sigw(i)^3+(9*pi*pi/8*(Hw/c/2)^2/Esp^2-3)*dg*(sigw(i))^2+3*dg^2*sigw(i)-dg^3);
        flag=i;
    end
end
sigw(flag)