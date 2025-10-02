function [CL,CH]=CLCH_Cal(kf,H,semi_angle,Ew,miu_w,Fg,dg) %kkk,scale factor
%kf is fracture toughness in unit of MPa m^(1/2)
%H is the hardness in Pa
%semi_angle is the semi angle of the indenter
%E is the young's modulus, miu is the poisson ratio
%%%%%% example:kf=4;H=2000*9.8*1e6;semi_angle=pi/2-pi/8;E=410e9;miu_w=0.28;Fg=2e-6;
C1=0.226; 
ss=cot(semi_angle); %for vicker
hi=pH_Cal(H,dg,Fg,miu_w,Ew);
a=sqrt(dg*hi); %contact radius; for sphere indenter, in mm
fai=asin(2*a/dg);
ss=pi/12*(3*tan(fai/2)+(tan(fai/2))^2); %for sphere indenter
%according to vicker indenter (experiment and theoretical invesigation...)
% CL=C1*(ss)^(5/6)*(Ew^(3/4)/H/(kf*1e6))^0.5*Fg^(5/8)*1e3;% m transfer to mm
% CH=C2*(ss)^(1/3)*Ew^(1/2)/H*Fg^(1/2)*1e3;% m transfer to mm
b=a*(Ew/H)^0.5*ss^(1/3);% equivalent to above, width of hemisphere plastic zone, in mm, ref:Elastic/Plastic Indentation Damage in Ceramics 
CH=b; %about 1.3um
% c=(0.032*sqrt(Ew/H)*(ss)^(2/3)*Fg/(kf*1e6))^(2/3)*1e3;%about 200 nm, %m transfer to mm, early radial crack? ref: the compelling case for..
% CL=c; 
% CL=4.5*a/2.718^(kf*1e6/(0.079*Fg/(a/1000)^1.5));% about 280nm,radial crack, ref:vickers indentation fracture toughness test
CL=4.5*a/2.718^(kf*1e6/(0.079/0.4636*Fg/(a/1000)^1.5));%about 1.2um if cosidering P=H*a^2 rather than Hv=0.4636*H*a^2,radial crack
% Cm=(4*0.064*Fg/(kf*1e6))^(2/3)*1000; % about 999nm,ref: a predictive model of the critical undeformed chip thickness
% CL=Cm/7; Vol=CL*CL*pi/2;
% b=sqrt(Fg/pi/H)*sqrt(3*(1-2*miu_w)/(5-4*miu_w)+2*sqrt(3)/pi/(5-4*miu_w)*Ew/(H/3)*ss) %stress field based, about 0.5um, ref: a new analytical model for estimation of scratch induced damage.
end