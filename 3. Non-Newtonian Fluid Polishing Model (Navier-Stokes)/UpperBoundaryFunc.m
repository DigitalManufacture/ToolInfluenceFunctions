function k=UpperBoundaryFunc(i,j,dx,dy,dz,lx,ly,lz,R,Rxy,gap,p_star,Ee)
    xbd=(i-1)*dx-lx/2;
    ybd=(j-1)*dy-ly/2;
        
    AA=0;Nx=floor(lx/dx);Ny=floor(ly/dy);Nz=floor(lz/dz);
        for ii=1:Nx
                for jj=1:Ny
                    xbd_ij=(ii-1)*dx-lx/2;
                    ybd_ij=(jj-1)*dy-ly/2;
                    rr=(xbd-xbd_ij)^2+(ybd-ybd_ij)^2;
                    if xbd_ij^2+ybd_ij^2<=Rxy^2
                       zbd_ij=R-sqrt(R^2-xbd_ij^2-ybd_ij^2)+gap;
                    else
                        zbd_ij=lz;
                    end
                    kk=round(zbd_ij/dz)+1;
                    if rr~=0
%                     if rr>=(deltax*2)^2
                       AA=AA+p_star(ii,jj,kk)*dx*dy/sqrt(rr)*2/pi/Ee;
                    end
                end
        end
        
        if xbd^2+ybd^2<=Rxy^2
            zbd=R-sqrt(R^2-xbd^2-ybd^2)+gap+AA;
        else
            zbd=lz;
        end
    k=round(zbd/dz)+1;
end