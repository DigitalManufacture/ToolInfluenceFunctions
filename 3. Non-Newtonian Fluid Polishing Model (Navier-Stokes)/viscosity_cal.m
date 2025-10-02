function visc=viscosity_cal(i,j,k,Us,Vs,Ws,dx,dy,dz,k0,n0)
    tmp=size(Us);
    maxsize=tmp(1)-1;
    if i==maxsize
        i=maxsize-1;
    end
    if i==0
        i=1;
    end
    tmp=size(Vs);
    maxsize=tmp(2)-1;
    if j==maxsize
        j=maxsize-1;
    end
    if j==0
        j=1;
    end
    tmp=size(Ws);
    maxsize=tmp(3)-1;
    if k==maxsize
        k=maxsize-1;
    end
    if k==0
        k=1;
    end
    
    GradV=[(Us(i+1,j,k)-Us(i,j,k))/dx (Vs(i+1,j,k)-Vs(i,j,k))/dx  (Ws(i+1,j,k)-Ws(i,j,k))/dx;
           (Us(i,j+1,k)-Us(i,j,k))/dy (Vs(i,j+1,k)-Vs(i,j,k))/dy  (Ws(i,j+1,k)-Ws(i,j,k))/dy;
           (Us(i,j,k+1)-Us(i,j,k))/dz (Vs(i,j,k+1)-Vs(i,j,k))/dz  (Ws(i,j,k+1)-Ws(i,j,k))/dz;
          ];
     kexi=(GradV+GradV')/2; 
     sum=0;
     for i=1:3
         for j=1:3
            sum=sum+kexi(i,j)^2;
         end
     end
     ShearRate=sqrt(sum);
     visc=k0*ShearRate^(n0-1)+0.06047; %from experimental fitting
     if visc==0
         visc=0.01;
     end
end