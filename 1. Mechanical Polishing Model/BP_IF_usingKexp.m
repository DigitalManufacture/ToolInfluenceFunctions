function [X Y Z r1]=BP_IF_usingKexp(Qt,wH,AlphaA,Kexp) %here wH in rpm, Alpha in degrees
% Qt=0.3; wH=800;AlphaA=20; %for test
Rt=11.1; % rubber radius
R2=12.1; % outside raius

miu_p=0.5; %poission ratio of polishing pad, polyuthan
% Qt=0.3; % tool offset
% wH=500*2*pi/60;% spindle speed
wH=wH*2*pi/60;% spindle speed

Ep=35e6;
r1=sqrt(Qt*R2);
ShoreH=60; % rubber shore hardness,90A for smooth plastic;
PoissonR=0.5; %rubber poisson ratio
Et=0.0981*(56+7.62336*ShoreH)/0.137505/(254-2.54*ShoreH)*1e6; %1e6:M
Eeq=1/((1-PoissonR^2)/Et);

Pmax=2*Eeq*sqrt(Qt)/pi/sqrt(Rt)/1e6;

% AlphaA=10*pi/180; %attack angle
AlphaA=AlphaA*pi/180; %attack angle
dg=1.5e-3; % grain size
fai=2/3; %percentage of wear flat
Ew=1.4e11; %ceramic Young's modulus, Pa
 
%for critical force cal
miu_w=0.28; %miu is the poisson ratio of ceramic
miu_grit=0.293; %diamond grit
Egrit=180e9; %1000GPa for diamond grit, 180GPa for CeO2

Ew=1/((1-miu_w^2)/Ew+(1-miu_grit^2)/Egrit); %Reudced Young's modulus, Pa

% sig_a=0.03; %planar radius of aspherity, mm
% A_a=pi*sig_a^2;
% R_a=(sig_a^2+l_a^2)/2/l_a;
R_a=0.1;
sig_a=0.035; %planar radius of aspherity, mm

delta_sig_a=2*sig_a;
Num_a=round(2*r1/delta_sig_a);
Dsum=Num_a^2/(2*r1)^2; %asperity number per area 1/mm^2


b1=pi*(3*R_a/4)^(2/3)*Dsum^(1/3);
Er=((1-miu_p^2)/Ep+(1-miu_w^2)/Ew)^(-1);% Pa

x=linspace(-r1,r1,Num_a);
y=linspace(-r1,r1,Num_a);
Pxy=NaN(length(x),length(x));

Pij=NaN(length(x),length(x));


MRR_th=NaN(length(x),length(x));
Vr=NaN(length(x),length(x));
bb=NaN(length(x),length(x));


for i=1:length(x)
    for j=1:length(y)
        if x(i)^2+y(j)^2<r1^2
           r_inspot=sqrt(x(i)^2+y(j)^2);
           Pxy(i,j)=Pr_Cal(Pmax,r1,r_inspot); %MPa
           bb(i,j)=b1*(Pxy(i,j)/(Er/1e6))^(2/3);
           Pij(i,j)=1/bb(i,j)*Pxy(i,j);%MPa
           %Pij(i,j)=Pxy(i,j);%MPa for smooth Pad
            
           Vr(i,j)=Vel_Cal(x(i),y(j),AlphaA,R2,Qt,wH);
           
          
        MRR_th(i,j)=Kexp*Pij(i,j)*Vr(i,j)*1e6;
        end
    end
end
[X,Y]=meshgrid(x,y);
Z=-MRR_th;
% [max(max(Pij)) trap_pressure p_cri]
% mesh(X,Y,Z)
for i=1:length(Z(:,1))
    for j=1:length(Z(:,1))
        if isnan(Z(i,j))
            Z(i,j)=0;
        end
    end
end
end

