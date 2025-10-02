%% need to calculate stress first
cc=30e-6; %g/mm^3, from 30g/L
dg_cluster=3.5e-3; %dg+1.5/2*2
dg=2.0e-3; %2 um grit size FO#6000
Ew=145e9; miu_w=0.4; %530 Alumina, 200 Nickel
H=123e6*9.8; %Yield strength 440MPa, refer to: Nickel as an Alternative Automotive Body Materials
% H=440e6*3;
Rou_g=3.65e-3; % g/mm^3
delta_g=(4/3*pi*(dg/2)^3*Rou_g/cc)^(1/3)*1;   %mm,uniform interval between grits
kk=dg_cluster/2/1000/dz; %correction coefficient due to the sparse grid

MRR=Us*0; MRR_h=Us*0; Pp=Us*0; h=Us*0;
for i=1:imax-1
    for j=1:jmax-1
        for k=1:kmax-1
            visc_all(i,j,k)=viscosity_cal(i,j,k,Us,Vs,Ws,dx,dy,dz,k0,n0);
            GradV=[(Us(i+1,j,k)-Us(i,j,k))/dx (Vs(i+1,j,k)-Vs(i,j,k))/dx  (Ws(i+1,j,k)-Ws(i,j,k))/dx;
                   (Us(i,j+1,k)-Us(i,j,k))/dy (Vs(i,j+1,k)-Vs(i,j,k))/dy  (Ws(i,j+1,k)-Ws(i,j,k))/dy;
                   (Us(i,j,k+1)-Us(i,j,k))/dz (Vs(i,j,k+1)-Vs(i,j,k))/dz  (Ws(i,j,k+1)-Ws(i,j,k))/dz;
                  ];
            kexi=(GradV+GradV')/2; 
            stressTensor=[-Ps(i,j,k) 0 0;0 -Ps(i,j,k) 0; 0 0 -Ps(i,j,k)]+visc_all(i,j,k)*kexi*2;
            
            Pp(i,j,k)=sqrt(1/2)*sqrt((stressTensor(1,1)-stressTensor(2,2))^2+(stressTensor(2,2)-stressTensor(3,3))^2+(stressTensor(1,1)-stressTensor(3,3))^2+6*(stressTensor(1,2)^2+stressTensor(2,3)^2+stressTensor(3,1)^2));
            Fg=Pp(i,j,k)*pi*(dg_cluster*1e-3)^2;
            h(i,j,k)=pH_Cal(H,dg,Fg,miu_w,Ew); %mm, Ew unit is Pa
            if h(i,j,k)<0; %remove infiinite small negative value
                h(i,j,k)=0;
            end
            MRR(i,j,k)=4/3*sqrt(2)*(dg/2)^0.5*sqrt(h(i,j,k)^3)*sqrt(Us(i,j,k)^2+Vs(i,j,k)^2)*1000*60; %mm^3/min
            MRR_h(i,j,k)=MRR(i,j,k)/delta_g/delta_g;
        end
    end
end
xp=linspace(-lx/2,lx/2,imax+1); yp=linspace(-ly/2,ly/2,jmax);
[X,Y]=meshgrid(xp,yp);   %vector plot

mesh(X'*1000,Y'*1000,-MRR_h(:,:,2)*kk);
figure;mesh(X'*1000,Y'*1000,h(:,:,2)); %depth in mm
figure;mesh(X'*1000,Y'*1000,Pp(:,:,2)*pi*(dg_cluster*1e-3)^2); % force on grit in N