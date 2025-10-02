function [u,v,w] = updateVelocity(imax,jmax,kmax,dx,dy,dz,lx,ly,lz,R,Rxy,gap,Ag,wH,u_star,v_star,w_star,p_prime,p,Ee,d_u,d_v,d_w,k0,n0,k_top1)
w = zeros(imax,jmax,kmax+1);
v = zeros(imax,jmax+1,kmax);
u = zeros(imax+1,jmax,kmax);

%update interior nodes of u and v
for i=2:imax
    for j=2:jmax-1
        for k=2:kmax-1
            u(i,j,k) = u_star(i,j,k) + d_u(i,j,k)*(p_prime(i-1,j,k)-p_prime(i,j,k));
        end
    end
end

for i=2:imax-1
    for j=2:jmax
        for k=2:kmax-1
            v(i,j,k) = v_star(i,j,k) + d_v(i,j,k)*(p_prime(i,j-1,k)-p_prime(i,j,k));
        end
    end
end

for i=2:imax-1
    for j=2:jmax-1
        for k=2:kmax
            w(i,j,k) = w_star(i,j,k) + d_w(i,j,k)*(p_prime(i,j,k-1)-p_prime(i,j,k));
        end
    end
end

%update BCs
% u(1,1:jmax,1:kmax) = -u(2,1:jmax,1:kmax); %left wall
% u(imax+1,1:jmax,1:kmax) = -u(imax,1:jmax,1:kmax); %right wall
% u(1:imax+1,jmax,1:kmax) = 0.0; %north wall
% u(1:imax+1,1,1:kmax) = 0.0; %south wall

u(1:imax+1,1:jmax, 1) = 0.0; %bottom wall
for i=1:imax+1 %top wall 
   for j=1:jmax
       %k_top=UpperBoundaryFunc(i,j,dx,dy,dz,lx,ly,lz,R,Rxy,gap,p,Ee);
       %k_top=UpperBoundaryFunc_Pp(i,j,dx,dy,dz,lx,ly,lz,R,Rxy,gap,p,Ee,u,v,w,k0,n0); %Principal stress induce deforma
       k_top=k_top1(i,j);
       [Umax, Vmax, Wmax]=UVW_boudary(i,j,lx,ly,R,Ag,wH,dx,dy);
        if k_top<kmax
            u(i,j,k_top) = Umax; 
        end
   end
end


% v(1:imax,1,1:kmax) = -v(1:imax,2,1:kmax); %south wall
% v(1:imax,jmax+1,1:kmax) = -v(1:imax,jmax,1:kmax); %north wall
% v(1,1:jmax+1,1:kmax) = 0.0; %left wall
% v(imax,1:jmax+1,1:kmax) = 0.0; %right wall

v(1:imax,1:jmax+1, 1) = 0.0; %bottom wall
for i=1:imax %top wall 
   for j=1:jmax+1
       %k_top=UpperBoundaryFunc(i,j,dx,dy,dz,lx,ly,lz,R,Rxy,gap,p,Ee);
       %k_top=UpperBoundaryFunc_Pp(i,j,dx,dy,dz,lx,ly,lz,R,Rxy,gap,p,Ee,u,v,w,k0,n0); %Principal stress induce deforma
       k_top=k_top1(i,j);
       [Umax, Vmax, Wmax]=UVW_boudary(i,j,lx,ly,R,Ag,wH,dx,dy);
       if k_top<kmax
          v(i,j,k_top) = Vmax; 
       end
   end
end


% w(1:imax,1,1:kmax) = 0; %south wall
% w(1:imax,jmax,1:kmax) = 0; %north wall
% w(1,1:jmax,1:kmax) = 0.0; %left wall
% w(imax,1:jmax,1:kmax) = 0.0; %right wall

w(1:imax,1:jmax, 1) = 0.0; %bottom wall
for i=1:imax %top wall 
   for j=1:jmax
       %k_top=UpperBoundaryFunc(i,j,dx,dy,dz,lx,ly,lz,R,Rxy,gap,p,Ee);
       %k_top=UpperBoundaryFunc_Pp(i,j,dx,dy,dz,lx,ly,lz,R,Rxy,gap,p,Ee,u,v,w,k0,n0); %Principal stress induce deforma
       k_top=k_top1(i,j);
       [Umax, Vmax, Wmax]=UVW_boudary(i,j,lx,ly,R,Ag,wH,dx,dy);
       if k_top<kmax
          w(i,j,k_top) = Wmax; 
       end
   end
end

return
end
