function [X Y Z r1]=BP_IF(Qt,wH,AlphaA) %here wH in rpm, Alpha in degrees
Rt=11.1; % rubber radius
R2=12.1; % outside raius
Th=1; %thickness of Polishing pad
Ep=3.5e9; %young's modulus of polishing pad, polyuthan
miu_p=0.5; %poission ratio of polishing pad, polyuthan
% Qt=0.3; % tool offset
% wH=500*2*pi/60;% spindle speed
wH=wH*2*pi/60;% spindle speed

r0=sqrt(Qt*Rt); %Hertz contact
r1=sqrt(Qt*R2);
ShoreH=60; % rubber shore hardness
PoissonR=0.5; %rubber poisson ratio
Et=0.0981*(56+7.62336*ShoreH)/0.137505/(254-2.54*ShoreH)*1e6; %1e6:M
Eeq=1/((1-PoissonR^2)/Et);
F1=4/3*Eeq/1e6*sqrt(Qt/Rt)*(r0)^2; 
c0=1/90; c1=1/90; %from experiments
F2=(c0*(Qt/Th)^0.5+c1*(Qt/Th)^2)*Ep/1e6*Th^3/R2;
Ftotal=F1+F2;
Eeff=3*Ftotal/4/sqrt(R2)/Qt^1.5; %MPa
Pmax=2*Eeff*sqrt(Qt)/pi/sqrt(R2);

% AlphaA=10*pi/180; %attack angle
AlphaA=AlphaA*pi/180; %attack angle
dg=1.5e-3; % grain size
fai=2/3; %percentage of wear flat
Ew=4.1e11; %ceramic Young's modulus, Pa
 
%for critical force cal
alpha=1/4; % dependent on indenter geometry, Brinell 1/4, vicker 2/pi
nt=1;theta=0.2; %coeffients
kf=4; %fracture toughness, MPa m^(1/2)
kkk=0.05; %scale factor
sc=0.6; %soften coefficient due to the CMP using CeO2
H=2000*9.8*1e6*sc; %viker hardness in Pa
sigmay=300*sc;  %yield stress in MPa
miu_w=0.28; %miu is the poisson ratio of ceramic
semi_angle=pi/2-pi/8; %semi_angle is the semi angle of the indenter
miu_plow=2/pi*(1/fai^2*asin(fai)-(1/fai^2-1)^0.5); %plowing coefficient

%define the pad asperity
l_a=0.02; %asperity height,mm, which is equal to the PV value of pad roughness
rou_percent=0.005; %0.005 of the total area are the asperities
delta_sig_a=0.03;
A_a=rou_percent*pi*r1^2/(2*r1/delta_sig_a)^2; %area per aspherity
Num_a=round(2*r1/delta_sig_a);
delta1=5e-6;% approximately roughness (base, assume proportion to pressure)

rou_g=3.9*1e-3; %g/mm^3
w_density=40; %g/L
G=w_density*1e-6/(rou_g*4/3*pi*(dg/2)^3);%concentration, number per mm^3

Er=((1-miu_p^2)/Ep+(1-miu_w^2)/Ew)^(-1);% Pa

x=linspace(-r1,r1,Num_a);
y=linspace(-r1,r1,Num_a);
Pxy=NaN(length(x),length(x));
Fgxy=NaN(length(x),length(x));
Pij=NaN(length(x),length(x));
N=NaN(length(x),length(x));
ch=NaN(length(x),length(x)); cl=NaN(length(x),length(x));
As=NaN(length(x),length(x));
MRR=NaN(length(x),length(x)); MRR_th=NaN(length(x),length(x));
da=NaN(length(x),length(x));
Vr=NaN(length(x),length(x));
bb=NaN(length(x),length(x));

[F_cri,h_cri]=Critical_Cal(alpha,nt,theta,kf,H,kkk,sigmay,dg);

for i=1:length(x)
    for j=1:length(y)
        if x(i)^2+y(j)^2<r1^2
           r_inspot=sqrt(x(i)^2+y(j)^2);
           Pxy(i,j)=Pr_Cal(Pmax,r1,r_inspot); %MPa
           bb(i,j)=0.25*Pxy(i,j)*dg/(H/(kf*1e6))/(delta1*Pxy(i,j)/Pmax); %assume delta1(roughness) is proportinonal to the pressure
           Pij(i,j)=1/bb(i,j)*Pxy(i,j);%MPa
           N(i,j)=G*bb(i,j)*l_a*A_a;
           Fgxy(i,j)=Pij(i,j)*pi*(dg/2)^2;
           if Fgxy(i,j)<=F_cri
               h(i,j)=pH_Cal(sigmay,dg,Fgxy(i,j));
               As(i,j)=4/3*sqrt(2)*sqrt(dg/2)*(h(i,j))^(3/2);
           else               
               [CL,CH,CS]=CLCH_Cal(kf,H,semi_angle,Ew,miu_w,Fgxy(i,j),dg,sigmay);
               ch(i,j)=CH; cl(i,j)=CL;
               As(i,j)=2*CL*CH;
           end
           Vr(i,j)=Vel_Cal(x(i),y(j),AlphaA,R2,Qt,wH);
           MRR(i,j)=N(i,j)*As(i,j)*Vr(i,j);
           MRR_th(i,j)=MRR(i,j)/pi/(bb(i,j)*A_a);
           MRR_th(i,j)=MRR(i,j)/pi/(A_a);
%            MRR_th(i,j)=MRR(i,j)/(x(2)-x(1))^2;
        end
    end
end
[X,Y]=meshgrid(x,y);
Z=-MRR_th;
end
