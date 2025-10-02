function [div]=checkDivergenceFree(imax,jmax,kmax,dx,dy,dz,u,v,w)

div=zeros(imax,jmax,kmax);

for i=1:imax
    for j=1:jmax
        for k=1:kmax
            div(i,j) = (1/dx)*(u(i,j,k)-u(i+1,j,k)) + (1/dy)*(v(i,j,k)-v(i,j+1,k))+ (1/dz)*(w(i,j,k)-w(i,j,k+1)); 
        end
    end
end

return
end