function [bp] = get_rhs(imax,jmax,kmax,dx,dy,dz,rho,u_star,v_star,w_star)

% RHS is the same for all nodes except the p_prime(1,1)
% because p(1,1) is set to be zero, it has no pressure correction

Ae=dy*dz; An=dx*dz; At=dx*dy;
Aw=Ae; As=An; Ab=At;
bp=zeros(imax,jmax,kmax);

for i=1:imax
    for j=1:jmax
        for k=1:kmax
            bp(i,j,k) = rho * (u_star(i,j,k)*Ae - u_star(i+1,j,k)*Ae + v_star(i,j,k)*An - v_star(i,j+1,k)*An+ w_star(i,j,k)*At - w_star(i,j,k+1)*At); 
        end
    end
end

% modify for p_prime(1,1)
bp(1,1,1) = 0;

return 
end
