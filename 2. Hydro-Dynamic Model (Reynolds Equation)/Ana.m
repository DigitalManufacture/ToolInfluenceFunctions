AlphaA=30*pi/180; %Precess angle
R2=22.0e-3; %tool radius
Qt=0.15e-3; % offset
nt=0.002; %viscosity
wH=12000*2*pi/60; %spindle speed
Ee=17e6; %33MPa  Tool Young's modulus
th=1.5e-3; %steady water thickness, m
Redge=sqrt(2*R2*th-th^2); %boundary
N=100; %res
[X Y]=meshgrid(linspace(-Redge,Redge,N),linspace(-Redge,Redge,N));
[U V]=Vel_Cal_UV(X,Y,AlphaA,R2,Qt,wH); %mm/s
V=abs(V);
% U=mean(mean(U))*ones(size(U)); %for test
% V=0*U;

delta=Redge*2/N;

con_a0=sqrt(R2*Qt);
PXY0=X;
for i=1:N
    for j=1:N
        if (X(i,j)^2+Y(i,j)^2<=con_a0^2)
            PXY0(i,j)=4/3*Ee*R2^0.5*Qt^1.5*3/2/pi/con_a0^2*(1-((X(i,j)^2+Y(i,j)^2))/con_a0^2);
        else
            PXY0(i,j)=0;
        end
    end
end

HXY_init=0*X;
H00=0.000001e-3;
kk=1.5;

for i=1:N
    for j=1:N
        if X(i,j)^2+Y(i,j)^2>con_a0^2
           HXY_init(i,j)=H00+R2-sqrt(R2^2-X(i,j)^2-Y(i,j)^2)-(R2-sqrt(R2^2-con_a0^2));
        else
            HXY_init(i,j)=H00;
        end 
     end
end
    
ERR=0.0002;
PK=PXY0*0;
ki=1;
while ki>0
    P=PK;
%     HXY0=hxy(N,R2,Ee,X,Y,con_a0,P,Redge);
HXY0=hxy_update(N,Ee,X,Y,P,Redge,HXY_init);
A=HXY0.^3; B=A;
C=3*HXY0.^2.*DIFF(HXY0,1,delta);
D=3*HXY0.^2.*DIFF(HXY0,2,delta);
E=6*U.*DIFF(HXY0,1,delta)*nt+6*V.*DIFF(HXY0,2,delta)*nt;

K=2*(A/delta/delta+B/delta/delta);
Cn=(B/delta/delta+D/2/delta)./K;
Cs=(B/delta/delta-D/2/delta)./K;
Ce=(A/delta/delta+C/2/delta)./K;
Cw=(A/delta/delta-C/2/delta)./K;
G=-E./K;

for i=2:N-1
    for j=2:N-1

        PK(i,j)=Cn(i,j)*P(i,j+1)+Cs(i,j)*P(i,j-1)+Ce(i,j)*P(i+1,j)+Cw(i,j)*P(i-1,j)+G(i,j);

        if PK(i,j)<0
            PK(i,j)=0;
            break;
        end
    end
end

sumP=sum(sum(P));
errP=sum(sum(P-PK));

ki=ki+1;
Fd=sum(sum(PK*delta*delta));
[ki abs(errP/sumP) Fd]

if abs(errP/sumP)<ERR
    break;
end
end

%%
HXY_20k=HXY0;
for i=1:N
    for j=1:N
        if X(i,j)^2+Y(i,j)^2>(Redge)^2
            HXY_20k(i,j)=NaN;
        end
    end
end
figure;mesh(X+0e-3,Y,HXY_20k)
set(gcf,'color','w')

