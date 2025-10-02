function [trap_depth trap_pressure]=trap_depth_Cal(H,dg,Ep,Hp,Cal_k) %kkk,scale factor
c=3;% relationship between hardness and yield stress
sigw=linspace(0,dg/2,1000); %estimation to reduce calculation
minA=abs(-dg^3);
flag=1;
pellet_elastic=(3*pi*0.4*Hp/4/Ep)^2*(dg/2);%refer eq. 11, a micro-contact and wear model...
% if pellet_elastic>dg
tg_s=(3*pi*0.4*Hp/4/Ep)^2*(dg/2);

for i=1:length(sigw)
    if abs(sigw(i)^3+(9*pi*pi/8*(H/c/2*Cal_k)^2/Ep^2-3)*dg*(sigw(i))^2+3*dg^2*sigw(i)-dg^3)<minA
        minA=abs(sigw(i)^3+(9*pi*pi/8*(H/c/2*Cal_k)^2/Ep^2-3)*dg*(sigw(i))^2+3*dg^2*sigw(i)-dg^3);
        flag=i;
    end
end
% end
if flag~=1
trap_depth=sigw(flag);
if dg-trap_depth<=tg_s
trap_pressure=H/1e6/c*pi*(dg/2)*trap_depth*Cal_k/pi/dg/dg*4; %Mpa
 [1 tg_s];  %for test
else
    [0 tg_s]; %for test
    trap_pressure=4*Hp*Cal_k*H/(2*Hp*c+Cal_k*H)/1e6; %Mpa
end

end
end