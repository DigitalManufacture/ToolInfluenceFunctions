function [F_cri,p_cri,h_cri,p_cri_elastic,h_cri_elastic]=Critical_Cal(alpha,nt,theta,kf,H,dg,miu_w,Ew,miu_grit,Egrit,tensile) %kkk,scale factor
%kf is fracture toughness in unit of MPa m^(1/2)
%H is the hardness in Pa
%alpha is dependent on indenter geometry, Brinell 1/4, vicker 2/pi
%sigmay is the yield strength in unit of MPa
%%%%%% example:alpha=1/4;nt=1;theta=0.2;kf=4;H=2000*9.8*1e6;,kkk=0.05;sigmay=300;dg=0.0015;
%ref:prediction and validation...(AMT), seems not reasonable
% kkk=1/2; %scale factor
% F_cri=kkk*54.47*alpha/nt^2/theta^4*(kf*1e6/H)^3*kf*1e6;  %about 0.017N
% h_cri=pH_Cal(H,dg,F_cri,miu_w,Ew);
%ref:spherical indentation of porous ceramics, seems not applicable
% Tao=6.5;%fracture energy, J/m^2
% k=9/16*(1-miu_w^2+(1-miu_grit^2)*Ew/Egrit);
% F_cri=1.17e5*Tao*k*dg/2*1e-3; %depends on the grit size
% F_cri=0.00042; %for test
% h_cri=pH_Cal(H,dg,F_cri,miu_w,Ew);

%ref:deformation and fracture of rocks loaded with speherical indenter
c=3; %relationship betwen yield stress and hardness, elastic-plastic solid is 3
% c value also refer:MECHANICAL PROPERTIES, INDENTATION AND DYNAMIC YIELD STRESS OF CERAMIC TARGETS
dafai=(1.1*pi/c)^3*(3*(1-miu_w)^2/4)^2;
F_cri_elastic=dafai*H^3/Ew^2*(dg/2*1e-3)^2; %in N, threshold to start yield
a=(3*(1-miu_w^2)*F_cri_elastic*dg/2*1e-3/4/Ew)^(1/3)*1e3; %in mm, contact radius
h_cri_elastic=a^2/(dg/2);

% kk=0.12;  % if kk=0.12 and tensile=8000e6 could match, equivalent to:
% Pmax_grain=F/pi/a^2>Sigma_t  ??
kk=2.75;
F_cri=(kk*tensile*2*pi/(1-2*miu_w))^3*(3*(1-miu_w^2)*dg*1e-3/8/Ew)^2;
% F_cri=F_cri/2;
p_cri=F_cri/pi/dg/dg*4; %Mpa crtical pressure on the grit to generate crack
h_cri=F_cri/(H/1e6/c)/pi/(dg/2)/Cal_k_coe(H,dg,miu_w,Ew);% based on plasiticiy, H/c is the yield stress
p_cri_elastic=F_cri_elastic/2/pi/dg/dg*4;
end