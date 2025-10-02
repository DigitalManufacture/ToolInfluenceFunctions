Qt=0.4; wH=800;AlphaA=30; %for test
Rt=11.1; % rubber radius
R2=12.1; % outside raius
Th=1; %thickness of Polishing pad
Ep=35e6; %500e6*5 resin for test, young's modulus of polishing pad, polyuthan, %%%%%%affect bb
Hp=90e6; %pellet hardness
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
c0=1/20; c1=1/10; %from experiments
F2=(c0*(Qt/Th)^0.5+c1*(Qt/Th)^2)*Ep/1e6*Th^3/R2;
Ftotal=F1+F2;
Eeff=3*Ftotal/4/sqrt(R2)/Qt^1.5; %MPa from exp
Eeff=10 ; %MPa
Pmax=2*Eeff*sqrt(Qt)/pi/sqrt(R2);

% AlphaA=10*pi/180; %attack angle
AlphaA=AlphaA*pi/180; %attack angle
dg=9e-3; % grain size
fai=2/3; %percentage of wear flat
Ew=1.5e11; %ceramic Young's modulus, Pa
 
%for critical force cal
alpha=1; % dependent on indenter geometry, Brinell 1, vicker 2/pi
nt=1;theta=0.2; %coeffients
kf=1.13; %fracture toughness, MPa m^(1/2)
sc=0.6; %soften coefficient due to the CMP using CeO2
H=1100*9.8*1e6; %viker hardness in Pa
c=3; %yield stress=H/c
tensile=1900e6; %workpiece tensile stress in PA
miu_w=0.28; %miu is the poisson ratio of ceramic
semi_angle=pi/2-pi/8; %semi_angle is the semi angle of the indenter
miu_plow=2/pi*(1/fai^2*asin(fai)-(1/fai^2-1)^0.5); %plowing coefficient
miu_grit=0.07; %diamond grit
Egrit=1220e9; %1000GPa for diamond grit

%define the pad asperity
l_a=0.02; %asperity height,mm
% sig_a=0.03; %planar radius of aspherity, mm
% A_a=pi*sig_a^2;
% R_a=(sig_a^2+l_a^2)/2/l_a;
R_a=0.1;
sig_a=0.03; %planar radius of aspherity, mm
l_a=R_a-sqrt(R_a^2-sig_a^2);
A_a=pi*sig_a^2;

delta_sig_a=2*sig_a;
Num_a=round(2*r1/delta_sig_a);
Dsum=Num_a^2/(2*r1)^2; %asperity number per area 1/mm^2

rou_g=3.9*1e-3; %g/mm^3
w_density=40; %g/L
G=w_density*1e-6/(rou_g*4/3*pi*(dg/2)^3);%concentration, number per mm^3

b1=pi*(3*R_a/4)^(2/3)*Dsum^(1/3);
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
h=NaN(length(x),length(x));

[F_cri,p_cri,h_cri,F_cri_elastic,h_cri_elastic]=Critical_Cal(alpha,nt,theta,kf,H,dg,miu_w,Ew,miu_grit,Egrit,tensile);
[trap_depth trap_pressure]=trap_depth_Cal(H,dg,Ep,Hp);

for i=1:length(x)
    for j=1:length(y)
        if x(i)^2+y(j)^2<r1^2
           r_inspot=sqrt(x(i)^2+y(j)^2);
           Pxy(i,j)=Pr_Cal(Pmax,r1,r_inspot); %MPa
           bb(i,j)=b1*(Pxy(i,j)/(Er/1e6))^(2/3);
           Pij(i,j)=3* Pxy(i,j);%MPa
%            N(i,j)=G*(pi*l_a*(3*sig_a^2+l_a^2)/6);
           N(i,j)=G*bb(i,j)*A_a*l_a;
           N(i,j)=440;
            
           Fgxy(i,j)=Pij(i,j)*pi*(0.25)^2/N(i,j); 
           if Fgxy(i,j)<=F_cri
               h(i,j)=pH_Cal(H,dg,Fgxy(i,j),miu_w,Ew);
               As(i,j)=4/3*sqrt(2)*sqrt(dg/2)*(h(i,j))^(3/2);
           else
               1
               [CL,CH]=CLCH_Cal(kf,H,semi_angle,Ew,miu_w,Fgxy(i,j),dg);
               ch(i,j)=CH; cl(i,j)=CL;
               As(i,j)=2*CL*CH;
           end
           Vr(i,j)=Vel_Cal(x(i),y(j),AlphaA,R2,Qt,wH);
           MRR(i,j)=N(i,j)*As(i,j)*Vr(i,j)*4;
           da(i,j)=Pxy(i,j)*pi*sig_a^2/Pij(i,j)/pi/R_a;
           MRR_th(i,j)=MRR(i,j)/pi/r0/r0;
%            MRR_th(i,j)=MRR(i,j)/(x(2)-x(1))^2;
        end
    end
end
[X,Y]=meshgrid(x,y);
Z=-MRR_th;
[max(max(Pij)) trap_pressure p_cri]

