%% initialization
clc
clear all
format longG
% GRID SIZE AND OTHER PARAMETERS
%i runs along x-direction and j runs along y-direction 

wH=3000*2*pi/60; %spindle speed
Ag=20/180*pi; % precess angle
R=20e-3; %tool radius
gap=0.12e-3; 
depth=0.1e-3; %immersed depth
n0=3.224; k0=2.5e-9; % fitted from experiments
Ee=33e6; %Tool youngs modulus

imax=33;                        %grid size in x-direction 
jmax=33;                        %grid size in y-direction 
kmax=38;                        %grid size in z-direction 

max_iteration=6000; 
maxRes = 1000;
iteration = 1;
mu = 0.01;                     %static viscosity
rho = 1e3;                     %density
lz=gap+depth; Rxy=sqrt(R^2-(R-depth)^2);
lx=10e-3; ly=10e-3; 
dx=lx/(imax-1);					%dx,dy cell sizes along x and y directions
dy=ly/(jmax-1); 
dz=lz/(kmax-1); 
x=dx/2:dx:lx-dx/2;              %plot use
y=dy/2:dy:ly-dy/2; 
z=dz/2:dz:lz-dz/2; 
alphaP = 0.2;                   %pressure under-relaxation, 0.4 previous
alphaU = 0.2;                   %velocity under-relaxation, 0.4 previous
tol = 10e-4;

%   u_star, v_star are Intermediate velocities
%   u and v = Final velocities

%Variable declaration
p   = zeros(imax,jmax,kmax);             %   p = Pressure
p_star   = zeros(imax,jmax,kmax);        
p_prime = zeros(imax,jmax,kmax);         %   pressure correction 
rhsp = zeros(imax,jmax,kmax);            %   Right hand side vector of pressure correction equation
divergence = zeros(imax,jmax,kmax); 
k_top=ones(imax+1,jmax+1)*kmax; 
AA=zeros(imax+1,jmax+1); 

% Horizontal Velocity -----------
u_star = zeros(imax+1,jmax,kmax);
uold   = zeros(imax+1,jmax,kmax);
uRes   = zeros(imax+1,jmax,kmax);
u      = zeros(imax+1,jmax,kmax);
d_u    = zeros(imax+1,jmax,kmax);  %velocity orrection coefficient

%Vertical velocity
v_star = zeros(imax,jmax+1,kmax);
vold   = zeros(imax,jmax+1,kmax);
vRes   = zeros(imax,jmax+1,kmax);
v      = zeros(imax,jmax+1,kmax);
d_v    = zeros(imax,jmax+1,kmax);    %velocity orrection coefficient

%Vertical velocity
w_star = zeros(imax,jmax,kmax+1);
wold   = zeros(imax,jmax,kmax+1);
wRes   = zeros(imax,jmax,kmax+1);
w      = zeros(imax,jmax,kmax+1);
d_w    = zeros(imax,jmax,kmax+1);    %velocity orrection coefficient

%% ---------- iterations -------------------
Converge_maxRes=[];
maxRes=tol+1; % to start the iteration
while ( (iteration <= max_iteration) && (maxRes > tol) ) 
    iteration = iteration + 1;
    NS=10; TS=50;
    if iteration<NS
        kkk=0;   %to converge, rigid boundary first (kkk=0), then after stable, kkk=1 is applied
    else if iteration<TS+NS
        kkk=(iteration-NS)/TS;
        else 
            kkk=1;
        end
    end
    kkk=0;  %remove this line and adjust variables NS and TS, if considering the structural coupling
    for i=1:imax
        for j=1:jmax
            [k_top(i,j) AA(i,j)]=UpperBoundaryFunc_Pp(i,j,dx,dy,dz,lx,ly,lz,R,Rxy,gap,p,Ee,u,v,w,k0,n0,kkk); %Principal stress induce deforma
        end
    end
    [u_star,d_u] = u_momentum(imax,jmax,kmax,lx,ly,lz,dx,dy,dz,rho,u,v,w,Ag,R,Rxy,gap,wH,k0,n0,p_star,Ee,alphaU,k_top);   %%Solve u-momentum equation for intermediate velocity u_star 
    [v_star,d_v] = v_momentum(imax,jmax,kmax,lx,ly,lz,dx,dy,dz,rho,u,v,w,Ag,R,Rxy,gap,wH,k0,n0,p_star,Ee,alphaU,k_top);   %%Solve v-momentum equation for intermediate velocity v_star 
    [w_star,d_w] = w_momentum(imax,jmax,kmax,lx,ly,lz,dx,dy,dz,rho,u,v,w,Ag,R,Rxy,gap,wH,k0,n0,p_star,Ee,alphaU,k_top);   %%Solve w-momentum equation for intermediate velocity w_star 
    uold = u;
    vold = v; 
    wold = w; 
    [rhsp] = get_rhs(imax,jmax,kmax,dx,dy,dz,rho,u_star,v_star,w_star);                                 %%Calculate rhs vector of the Pressure Poisson matrix 
    [p,Pp] = pres_correct_interior(imax,jmax,kmax,dx,dy,dz,lx,ly,lz,R,Rxy,gap,rhsp,rho,d_u,d_v,d_w,p_star,Ee,alphaP,u,v,w,k0,n0,k_top);   %%Solve pressure correction implicitly and update pressure
    [u,v,w] = updateVelocity(imax,jmax,kmax,dx,dy,dz,lx,ly,lz,R,Rxy,gap,Ag,wH,u_star,v_star,w_star,Pp,p,Ee,d_u,d_v,d_w,k0,n0,k_top);            %%Update velocity based on pressure correction
    [divergence]=checkDivergenceFree(imax,jmax,kmax,dx,dy,dz,u,v,w);                               %%check if velocity field is divergence free
    [XT,YT,ZT] = ToolBoundaryCheck(imax,jmax,dx,dy,dz,lx,ly,lz,R,Rxy,gap,p,Ee,u,v,w,k0,n0,k_top);
    p_star = p;                                                                          %%use p as p_star for the next iteration

    %find maximum residual in the domain
    vRes = abs(v - vold);
    uRes = abs(u - uold);
    wRes = abs(w - wold);

    maxRes_u = max(max(max(uRes)));
    maxRes_v = max(max(max(vRes)));
    maxRes_w = max(max(max(wRes)));
    
    meanRes_u = mean(mean(mean(uRes)));
    meanRes_v = mean(mean(mean(vRes)));
    meanRes_w = mean(mean(mean(wRes)));

    maxRes = max([maxRes_u, maxRes_v, maxRes_w]);
    meanRes = max([meanRes_u, meanRes_v, meanRes_w]);
    
    Converge_maxRes=[Converge_maxRes; maxRes meanRes];                                                                        %%Check for convergence 
    disp(['Iter = ',int2str(iteration),'; Res = ',num2str(maxRes),'; CenterDeformed (mm)= ',num2str(AA(round(imax/2),round(jmax/2))*kkk*1000),'; kkk=',num2str(kkk)]);
    if (maxRes > 10)
        disp('not going to converge!');
        break;
    end
end
%% calculate the equivalent stress
Ps=p; Us=u; Vs=v; Ws=w; Pp=p;
visc_all=Ps*0;
ss_all_1=Ps*0; ss_all_2=Ps*0; ss_all_3=Ps*0;
nr=[0 0 1]; %which face
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
            stressface=nr*stressTensor;
            ss_all_1(i,j,k)=stressface(1);
            ss_all_2(i,j,k)=stressface(2);
            ss_all_3(i,j,k)=stressface(3);
        end
    end
end
%figure;mesh(-p(:,:,10)/2+sqrt((-p(:,:,10)/2).^2+ss_all_1(:,:,10).^2+ss_all_2(:,:,10).^2+ss_all_3(:,:,10).^2));
figure;mesh(Pp(:,:,2))
% cp=3;
% kk=cp:imax-cp;
% surf(visc_all(kk,kk,10))
% surf(Us(kk,kk,10))
