function [u_star,d_u] = u_momentum(imax,jmax,kmax,lx,ly,lz,dx,dy,dz,rho,u,v,w,Ag,R,Rxy,gap,wH,k0,n0,p,Ee,alpha,k_top1)

u_star=zeros(imax+1,jmax,kmax);
d_u=zeros(imax+1,jmax,kmax);

Ae=dy*dz; An=dx*dz; At=dx*dy;
Aw=Ae; As=An; Ab=At;

A = @(F,D)( max(0, (1-0.1 * abs(F/D))^5 ) );

%%compute u_star
for i = 2:imax
    for j = 2:jmax-1
       %k_top=UpperBoundaryFunc(i,j,dx,dy,dz,lx,ly,lz,R,Rxy,gap,p,Ee); %dynamic pressure only induce deforma
       %k_top=UpperBoundaryFunc_Pp(i,j,dx,dy,dz,lx,ly,lz,R,Rxy,gap,p,Ee,u,v,w,k0,n0); %Principal stress induce deforma
       k_top=k_top1(i,j);
       for k=2:k_top-1
        De=Ae/(dx/2)*viscosity_cal(i,j,k-1,u,v,w,dx,dy,dz,k0,n0);
        Dw=Aw/(dx/2)*viscosity_cal(i-1,j,k-1,u,v,w,dx,dy,dz,k0,n0);
        Dn=1/2*An/(dy/2)*viscosity_cal(i,j,k-1,u,v,w,dx,dy,dz,k0,n0);
        Ds=1/2*As/(dy/2)*viscosity_cal(i,j-1,k-1,u,v,w,dx,dy,dz,k0,n0);
        Dt=1/2*At/(dz/2)*viscosity_cal(i,j,k,u,v,w,dx,dy,dz,k0,n0);
        Db=1/2*Ab/(dz/2)*viscosity_cal(i,j,k-1,u,v,w,dx,dy,dz,k0,n0);
        
        Fe  = .5*rho*Ae*(u(i+1,j,k)+u(i,j,k));                                     
        Fw  = .5*rho*Ae*(u(i-1,j,k)+u(i,j,k)); 
        Fn  = .5*rho*An*(v(i,j+1,k)+v(i-1,j+1,k)); 
        Fs  = .5*rho*An*(v(i,j,k)+v(i-1,j,k));
        Ft  = .5*rho*At*(w(i,j,k+1)+w(i,j,k+1)); 
        Fb  = .5*rho*At*(w(i,j,k)+w(i,j,k));
        
        aE = De * A(Fe,De) + max(-Fe,0);
        aW = Dw * A(Fw,Dw) + max(Fw,0);
        aN = Dn * A(Fn,Dn) + max(-Fn,0);
        aS = Ds * A(Fs,Ds) + max(Fs,0);
        aT = Dt * A(Ft,Dt) + max(-Ft,0);
        aB = Db * A(Fb,Db) + max(Fb,0);
        
        aP = aE + aW + aN + aS + aT + aB + (Fe-Fw) + (Fn-Fs) + (Ft-Fb);
        
        pressure_term = (p(i-1,j,k)-p(i,j,k)) * Ae;        
       
        u_star(i,j,k) = alpha/aP * ( (aE*u(i+1,j,k)+aW*u(i-1,j,k)+aN*u(i,j+1,k)+aS*u(i,j-1,k))+aT*u(i,j,k+1)+aB*u(i,j,k-1) + pressure_term ) + (1-alpha)*u(i,j,k);
        
        d_u(i,j,k) = alpha * Ae / aP;   %refer to Versteeg CFD book
        end
    end
end

%%set d_u for top and bottom BCs
%%they will be later used by the pressure correction equation 
%%they should not be zero, or BCs of pressure correction will get messed up
k = 1; %bottom
for i=2:imax
   for j=2:jmax
    
        De=Ae/(dx/2)*viscosity_cal(i,j,k-1,u,v,w,dx,dy,dz,k0,n0); % 
        Dw=Aw/(dx/2)*viscosity_cal(i-1,j,k-1,u,v,w,dx,dy,dz,k0,n0);
        Dn=1/2*An/(dy/2)*viscosity_cal(i,j,k-1,u,v,w,dx,dy,dz,k0,n0);
        Ds=1/2*As/(dy/2)*viscosity_cal(i,j-1,k-1,u,v,w,dx,dy,dz,k0,n0);
        Dt=1/2*At/(dz/2)*viscosity_cal(i,j,k,u,v,w,dx,dy,dz,k0,n0);
        Db=1/2*Ab/(dz/2)*viscosity_cal(i,j,k-1,u,v,w,dx,dy,dz,k0,n0);
        
        Fe  = .5*rho*Ae*(u(i+1,j,k)+u(i,j,k));                                     
        Fw  = .5*rho*Ae*(u(i-1,j,k)+u(i,j,k)); 
        Fn  = .5*rho*An*(v(i,j+1,k)+v(i-1,j+1,k)); 
        Fs  = .5*rho*An*(v(i,j,k)+v(i-1,j,k));
        Ft  = .5*rho*At*(w(i,j,k+1)+w(i-1,j,k+1)); 
        Fb  = 0;
        
        aE = De * A(Fe,De) + max(-Fe,0);
        aW = Dw * A(Fw,Dw) + max(Fw,0);
        aN = Dn * A(Fn,Dn) + max(-Fn,0);
        aS = Ds * A(Fs,Ds) + max(Fs,0);
        aT = Dt * A(Ft,Dt) + max(-Ft,0);
        aB = 0;
        
        aP = aE + aW + aN + aS + aT + aB + (Fe-Fw) + (Fn-Fs) + (Ft-Fb);
        
        d_u(i,j,k) = alpha * Ae / aP;
    end
end

%top
for i=2:imax
   for j=2:jmax
       
        %k=UpperBoundaryFunc(i,j,dx,dy,dz,lx,ly,lz,R,Rxy,gap,p,Ee);
        %k=UpperBoundaryFunc_Pp(i,j,dx,dy,dz,lx,ly,lz,R,Rxy,gap,p,Ee,u,v,w,k0,n0); %Principal stress induce deforma
        k=k_top1(i,j);
        De=Ae/(dx/2)*viscosity_cal(i,j,k-1,u,v,w,dx,dy,dz,k0,n0); % 
        Dw=Aw/(dx/2)*viscosity_cal(i-1,j,k-1,u,v,w,dx,dy,dz,k0,n0);
        Dn=1/2*An/(dy/2)*viscosity_cal(i,j,k-1,u,v,w,dx,dy,dz,k0,n0);
        Ds=1/2*As/(dy/2)*viscosity_cal(i,j-1,k-1,u,v,w,dx,dy,dz,k0,n0);
        Dt=1/2*At/(dz/2)*viscosity_cal(i,j,k,u,v,w,dx,dy,dz,k0,n0);
        Db=1/2*Ab/(dz/2)*viscosity_cal(i,j,k-1,u,v,w,dx,dy,dz,k0,n0);
        
        Fe  = .5*rho*Ae*(u(i+1,j,k)+u(i,j,k));                                     
        Fw  = .5*rho*Ae*(u(i-1,j,k)+u(i,j,k)); 
        Fn  = .5*rho*An*(v(i,j+1,k)+v(i-1,j+1,k)); 
        Fs  = .5*rho*An*(v(i,j,k)+v(i-1,j,k));
        Ft  = 0; 
        Fb  = .5*rho*At*(w(i,j,k)+w(i-1,j,k));
        
        aE = De * A(Fe,De) + max(-Fe,0);
        aW = Dw * A(Fw,Dw) + max(Fw,0);
        aN = Dn * A(Fn,Dn) + max(-Fn,0);
        aS = Ds * A(Fs,Ds) + max(Fs,0);
        aT = 0;
        aB = Db * A(Fb,Db) + max(Fb,0);
        
        aP = aE + aW + aN + aS + aT + aB + (Fe-Fw) + (Fn-Fs) + (Ft-Fb);
        
        d_u(i,j,k) = alpha * Ae / aP;
   end
end

j = 1; %south
for i=2:imax
   for k=2:kmax
    
        De=Ae/(dx/2)*viscosity_cal(i,j,k-1,u,v,w,dx,dy,dz,k0,n0); % k 
        Dw=Aw/(dx/2)*viscosity_cal(i-1,j,k-1,u,v,w,dx,dy,dz,k0,n0);
        Dn=1/2*An/(dy/2)*viscosity_cal(i,j,k-1,u,v,w,dx,dy,dz,k0,n0);
        Ds=1/2*As/(dy/2)*viscosity_cal(i,j-1,k-1,u,v,w,dx,dy,dz,k0,n0);
        Dt=1/2*At/(dz/2)*viscosity_cal(i,j,k,u,v,w,dx,dy,dz,k0,n0);
        Db=1/2*Ab/(dz/2)*viscosity_cal(i,j,k-1,u,v,w,dx,dy,dz,k0,n0);
        
        Fe  = .5*rho*Ae*(u(i+1,j,k)+u(i,j,k));                                     
        Fw  = .5*rho*Ae*(u(i-1,j,k)+u(i,j,k)); 
        Fn  = .5*rho*An*(v(i,j+1,k)+v(i-1,j+1,k)); 
        Fs  = 0;
        Ft  = .5*rho*At*(w(i,j,k+1)+w(i-1,j,k+1)); 
        Fb  = .5*rho*At*(w(i,j,k)+w(i-1,j,k));
        
        aE = De * A(Fe,De) + max(-Fe,0);
        aW = Dw * A(Fw,Dw) + max(Fw,0);
        aN = Dn * A(Fn,Dn) + max(-Fn,0);
        aS = 0;
        aT = Dt * A(Ft,Dt) + max(-Ft,0);
        aB = Db * A(Fb,Db) + max(Fb,0);
        
        aP = aE + aW + aN + aS + aT + aB + (Fe-Fw) + (Fn-Fs) + (Ft-Fb);
        
        d_u(i,j,k) = alpha * Ae / aP;
    end
end

j = jmax; %north
for i=2:imax
   for k=2:kmax
    
        De=Ae/(dx/2)*viscosity_cal(i,j,k-1,u,v,w,dx,dy,dz,k0,n0); % k
        Dw=Aw/(dx/2)*viscosity_cal(i-1,j,k-1,u,v,w,dx,dy,dz,k0,n0);
        Dn=1/2*An/(dy/2)*viscosity_cal(i,j,k-1,u,v,w,dx,dy,dz,k0,n0);
        Ds=1/2*As/(dy/2)*viscosity_cal(i,j-1,k-1,u,v,w,dx,dy,dz,k0,n0);
        Dt=1/2*At/(dz/2)*viscosity_cal(i,j,k,u,v,w,dx,dy,dz,k0,n0);
        Db=1/2*Ab/(dz/2)*viscosity_cal(i,j,k-1,u,v,w,dx,dy,dz,k0,n0);
        
        Fe  = .5*rho*Ae*(u(i+1,j,k)+u(i,j,k));                                     
        Fw  = .5*rho*Ae*(u(i-1,j,k)+u(i,j,k)); 
        Fn  = 0; 
        Fs  = .5*rho*An*(v(i,j,k)+v(i-1,j,k));
        Ft  = .5*rho*At*(w(i,j,k+1)+w(i-1,j,k+1)); 
        Fb  = .5*rho*At*(w(i,j,k)+w(i-1,j,k));
        
        aE = De * A(Fe,De) + max(-Fe,0);
        aW = Dw * A(Fw,Dw) + max(Fw,0);
        aN = 0;
        aS = Ds * A(Fs,Ds) + max(Fs,0);
        aT = Dt * A(Ft,Dt) + max(-Ft,0);
        aB = Db * A(Fb,Db) + max(Fb,0);
        
        aP = aE + aW + aN + aS + aT + aB + (Fe-Fw) + (Fn-Fs) + (Ft-Fb);
        
        d_u(i,j,k) = alpha * Ae / aP;
    end
end

%%Apply BCs
% u_star(1,1:jmax,1:kmax) = -u_star(2,1:jmax,1:kmax); %left wall
% u_star(imax+1,1:jmax,1:kmax) = -u_star(imax,1:jmax,1:kmax); %right wall
% u_star(1:imax+1,jmax,1:kmax) = 0.0; %north wall
% u_star(1:imax+1,1,1:kmax) = 0.0; %south wall

u_star(1:imax+1,1:jmax, 1) = 0.0; %bottom wall
for i=1:imax+1 %top wall 
   for j=1:jmax
       %k_top=UpperBoundaryFunc(i,j,dx,dy,dz,lx,ly,lz,R,Rxy,gap,p,Ee);
       %k_top=UpperBoundaryFunc_Pp(i,j,dx,dy,dz,lx,ly,lz,R,Rxy,gap,p,Ee,u,v,w,k0,n0); %Principal stress induce deforma
       k_top=k_top1(i,j);
       [Umax, Vmax, Wmax]=UVW_boudary(i,j,lx,ly,R,Ag,wH,dx,dy);
        if k_top<kmax
           u_star(i,j,k_top) = Umax; 
        end
   end
end

return 
end