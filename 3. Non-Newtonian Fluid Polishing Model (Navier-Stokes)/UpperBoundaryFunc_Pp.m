function [k AA]=UpperBoundaryFunc_Pp(i,j,dx,dy,dz,lx,ly,lz,R,Rxy,gap,p_star,Ee,u,v,w,k0,n0,kkk)

    xbd=(i-1)*dx-lx/2;
    ybd=(j-1)*dy-ly/2;
        
    AA=0;Nx=floor(lx/dx);Ny=floor(ly/dy);Nz=floor(lz/dz);
        for ii=1:Nx-1
                for jj=1:Ny-1
                    xbd_ij=(ii-1)*dx-lx/2;
                    ybd_ij=(jj-1)*dy-ly/2;
                    rr=(xbd-xbd_ij)^2+(ybd-ybd_ij)^2;
                    if xbd_ij^2+ybd_ij^2<=Rxy^2
                       zbd_ij=R-sqrt(R^2-xbd_ij^2-ybd_ij^2)+gap;
                    else
                        zbd_ij=lz;
                    end
                  if xbd_ij^2+ybd_ij^2<=Rxy^2
                    kk=round(zbd_ij/dz)-5;  %modified,original kk=round(zbd_ij/dz)+1
                    if kk>=round(lz/dz)   % Principal stress induce def
                        kk=round(lz/dz);
                    end
                    if kk<1
                        kk=1;
                    end
                    Us=u;Vs=v;Ws=w;
                    visc_iijjkk=viscosity_cal(ii,jj,kk,Us,Vs,Ws,dx,dy,dz,k0,n0);
            GradV=[(Us(ii+1,jj,kk)-Us(ii,jj,kk))/dx (Vs(ii+1,jj,kk)-Vs(ii,jj,kk))/dx  (Ws(ii+1,jj,kk)-Ws(ii,jj,kk))/dx;
                   (Us(ii,jj+1,kk)-Us(ii,jj,kk))/dy (Vs(ii,jj+1,kk)-Vs(ii,jj,kk))/dy  (Ws(ii,jj+1,kk)-Ws(ii,jj,kk))/dy;
                   (Us(ii,jj,kk+1)-Us(ii,jj,kk))/dz (Vs(ii,jj,kk+1)-Vs(ii,jj,kk))/dz  (Ws(ii,jj,kk+1)-Ws(ii,jj,kk))/dz;
                  ];
            kexi=(GradV+GradV')/2; 
            stressTensor=[-p_star(ii,jj,kk) 0 0;0 -p_star(ii,jj,kk) 0; 0 0 -p_star(ii,jj,kk)]+visc_iijjkk*kexi*2;
            Pp_iijjkk=sqrt(1/2)*sqrt((stressTensor(1,1)-stressTensor(2,2))^2+(stressTensor(2,2)-stressTensor(3,3))^2+(stressTensor(1,1)-stressTensor(3,3))^2+6*(stressTensor(1,2)^2+stressTensor(2,3)^2+stressTensor(3,1)^2));
                    if rr~=0
                       AA=AA+Pp_iijjkk*dx*dy/sqrt(rr)*2/pi/Ee;
                    end
                  end
%                     if rr~=0
%                        AA=AA+p_star(ii,jj,kk)*dx*dy/sqrt(rr)*2/pi/Ee;
%                     end
                end
        end
        
        if xbd^2+ybd^2<=Rxy^2
            zbd=R-sqrt(R^2-xbd^2-ybd^2)+gap;
        else
            zbd=lz;
        end
        if zbd>lz
            zbd=lz;
        end
    k=round(zbd/dz)+1+round(AA*kkk/dz);
    if k>=Nz+1
        k=Nz+1;
    end
end