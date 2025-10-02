function [v_star,d_v] = v_momentum(imax,jmax,kmax,lx,ly,lz,dx,dy,dz,rho,u,v,w,Ag,R,Rxy,gap,wH,k0,n0,p,Ee,alpha,k_top1)

v_star=zeros(imax,jmax+1,kmax);
d_v=zeros(imax,jmax+1,kmax);

Ae=dy*dz; An=dx*dz; At=dx*dy;
Aw=Ae; As=An; Ab=At;

A = @(F,D)( max(0, (1-0.1 * abs(F/D))^5 ) );

%%compute v_star
for i = 2:imax-1
    for j = 2:jmax
       %k_top=UpperBoundaryFunc(i,j,dx,dy,dz,lx,ly,lz,R,Rxy,gap,p,Ee);
       %k_top=UpperBoundaryFunc_Pp(i,j,dx,dy,dz,lx,ly,lz,R,Rxy,gap,p,Ee,u,v,w,k0,n0); %Principal stress induce deforma
       k_top=k_top1(i,j);
       for k=2:k_top-1
        De=1/2*Ae/(dx/2)*viscosity_cal(i,j,k-1,u,v,w,dx,dy,dz,k0,n0);
        Dw=1/2*Aw/(dx/2)*viscosity_cal(i-1,j,k-1,u,v,w,dx,dy,dz,k0,n0);
        Dn=An/(dy/2)*viscosity_cal(i,j,k-1,u,v,w,dx,dy,dz,k0,n0);
        Ds=As/(dy/2)*viscosity_cal(i,j-1,k-1,u,v,w,dx,dy,dz,k0,n0);
        Dt=1/2*At/(dz/2)*viscosity_cal(i,j,k,u,v,w,dx,dy,dz,k0,n0);
        Db=1/2*Ab/(dz/2)*viscosity_cal(i,j,k-1,u,v,w,dx,dy,dz,k0,n0);
        
        Fe  = .5*rho*Ae*(u(i+1,j,k)+u(i+1,j-1,k));                                     
        Fw  = .5*rho*Ae*(u(i,j,k)+u(i,j-1,k)); 
        Fn  = .5*rho*An*(v(i,j,k)+v(i,j+1,k)); 
        Fs  = .5*rho*An*(v(i,j-1,k)+v(i,j,k));
        Ft  = .5*rho*At*(w(i,j,k+1)+w(i,j,k+1)); 
        Fb  = .5*rho*At*(w(i,j,k)+w(i,j,k));
                
        aE = De * A(Fe,De) + max(-Fe,0);
        aW = Dw * A(Fw,Dw) + max(Fw,0);
        aN = Dn * A(Fn,Dn) + max(-Fn,0);
        aS = Ds * A(Fs,Ds) + max(Fs,0);
        aT = Dt * A(Ft,Dt) + max(-Ft,0);
        aB = Db * A(Fb,Db) + max(Fb,0);
        
        aP = aE + aW + aN + aS + aT + aB + (Fe-Fw) + (Fn-Fs) + (Ft-Fb);
        
        pressure_term = (p(i,j-1,k)-p(i,j,k)) * An;        
       
        v_star(i,j,k) = alpha/aP * ( (aE*v(i+1,j,k)+aW*v(i-1,j,k)+aN*v(i,j+1,k)+aS*v(i,j-1,k))+aT*v(i,j,k+1)+aB*v(i,j,k-1) + pressure_term ) + (1-alpha)*v(i,j,k);
        
        d_v(i,j,k) = alpha * An / aP;   %refer to Versteeg CFD book
        end
    end
end

%%set d_u for top and bottom BCs
%%they will be later used by the pressure correction equation 
%%they should not be zero, or BCs of pressure correction will get messed up
k = 1; %bottom
for i=2:imax
   for j=2:jmax
    
        De=1/2*Ae/(dx/2)*viscosity_cal(i,j,k-1,u,v,w,dx,dy,dz,k0,n0); % 
        Dw=1/2*Aw/(dx/2)*viscosity_cal(i-1,j,k-1,u,v,w,dx,dy,dz,k0,n0);
        Dn=An/(dy/2)*viscosity_cal(i,j,k-1,u,v,w,dx,dy,dz,k0,n0);
        Ds=As/(dy/2)*viscosity_cal(i,j-1,k-1,u,v,w,dx,dy,dz,k0,n0);
        Dt=1/2*At/(dz/2)*viscosity_cal(i,j,k,u,v,w,dx,dy,dz,k0,n0);
        Db=1/2*Ab/(dz/2)*viscosity_cal(i,j,k-1,u,v,w,dx,dy,dz,k0,n0);
        
        Fe  = .5*rho*Ae*(u(i+1,j,k)+u(i+1,j-1,k));                                     
        Fw  = .5*rho*Ae*(u(i,j,k)+u(i,j-1,k)); 
        Fn  = .5*rho*An*(v(i,j,k)+v(i,j+1,k)); 
        Fs  = .5*rho*An*(v(i,j-1,k)+v(i,j,k));
        Ft  = .5*rho*At*(w(i,j-1,k+1)+w(i,j,k+1)); 
        Fb  = 0;
        
        aE = De * A(Fe,De) + max(-Fe,0);
        aW = Dw * A(Fw,Dw) + max(Fw,0);
        aN = Dn * A(Fn,Dn) + max(-Fn,0);
        aS = Ds * A(Fs,Ds) + max(Fs,0);
        aT = Dt * A(Ft,Dt) + max(-Ft,0);
        aB = 0;
        
        aP = aE + aW + aN + aS + aT + aB + (Fe-Fw) + (Fn-Fs) + (Ft-Fb);
        
        d_v(i,j,k) = alpha * An / aP;
    end
end

%top
for i=2:imax
   for j=2:jmax
       
        %k=UpperBoundaryFunc(i,j,dx,dy,dz,lx,ly,lz,R,Rxy,gap,p,Ee);
        %k=UpperBoundaryFunc_Pp(i,j,dx,dy,dz,lx,ly,lz,R,Rxy,gap,p,Ee,u,v,w,k0,n0); %Principal stress induce deforma
        k=k_top1(i,j);
        De=1/2*Ae/(dx/2)*viscosity_cal(i,j,k-1,u,v,w,dx,dy,dz,k0,n0); % 
        Dw=1/2*Aw/(dx/2)*viscosity_cal(i-1,j,k-1,u,v,w,dx,dy,dz,k0,n0);
        Dn=An/(dy/2)*viscosity_cal(i,j,k-1,u,v,w,dx,dy,dz,k0,n0);
        Ds=As/(dy/2)*viscosity_cal(i,j-1,k-1,u,v,w,dx,dy,dz,k0,n0);
        Dt=1/2*At/(dz/2)*viscosity_cal(i,j,k,u,v,w,dx,dy,dz,k0,n0);
        Db=1/2*Ab/(dz/2)*viscosity_cal(i,j,k-1,u,v,w,dx,dy,dz,k0,n0);
        
        Fe  = .5*rho*Ae*(u(i+1,j,k)+u(i+1,j-1,k));                                     
        Fw  = .5*rho*Ae*(u(i,j,k)+u(i,j-1,k)); 
        Fn  = .5*rho*An*(v(i,j,k)+v(i,j+1,k)); 
        Fs  = .5*rho*An*(v(i,j-1,k)+v(i,j,k));
        Ft  = 0; 
        Fb  = .5*rho*At*(w(i,j-1,k)+w(i,j,k));
        
        aE = De * A(Fe,De) + max(-Fe,0);
        aW = Dw * A(Fw,Dw) + max(Fw,0);
        aN = Dn * A(Fn,Dn) + max(-Fn,0);
        aS = Ds * A(Fs,Ds) + max(Fs,0);
        aT = 0;
        aB = Db * A(Fb,Db) + max(Fb,0);
        
        aP = aE + aW + aN + aS + aT + aB + (Fe-Fw) + (Fn-Fs) + (Ft-Fb);
        
        d_v(i,j,k) = alpha * An / aP;
   end
end

i=1; %left
for j = 2:jmax
       %k_top=UpperBoundaryFunc(i,j,dx,dy,dz,lx,ly,lz,R,Rxy,gap,p,Ee);
       %k_top=UpperBoundaryFunc_Pp(i,j,dx,dy,dz,lx,ly,lz,R,Rxy,gap,p,Ee,u,v,w,k0,n0); %Principal stress induce deforma
       k_top=k_top1(i,j);
       for k=2:k_top-1
        De=1/2*Ae/(dx/2)*viscosity_cal(i,j,k-1,u,v,w,dx,dy,dz,k0,n0);
        Dw=1/2*Aw/(dx/2)*viscosity_cal(i-1,j,k-1,u,v,w,dx,dy,dz,k0,n0);
        Dn=An/(dy/2)*viscosity_cal(i,j,k-1,u,v,w,dx,dy,dz,k0,n0);
        Ds=As/(dy/2)*viscosity_cal(i,j-1,k-1,u,v,w,dx,dy,dz,k0,n0);
        Dt=1/2*At/(dz/2)*viscosity_cal(i,j,k,u,v,w,dx,dy,dz,k0,n0);
        Db=1/2*Ab/(dz/2)*viscosity_cal(i,j,k-1,u,v,w,dx,dy,dz,k0,n0);
        
        Fe  = .5*rho*Ae*(u(i+1,j,k)+u(i+1,j-1,k));                                     
        Fw  = 0; 
        Fn  = .5*rho*An*(v(i,j,k)+v(i,j+1,k)); 
        Fs  = .5*rho*An*(v(i,j-1,k)+v(i,j,k));
        Ft  = .5*rho*At*(w(i,j-1,k+1)+w(i,j,k+1)); 
        Fb  = .5*rho*At*(w(i,j-1,k)+w(i,j,k));
        
        aE = De * A(Fe,De) + max(-Fe,0);
        aW = 0;
        aN = Dn * A(Fn,Dn) + max(-Fn,0);
        aS = Ds * A(Fs,Ds) + max(Fs,0);
        aT = Dt * A(Ft,Dt) + max(-Ft,0);
        aB = Db * A(Fb,Db) + max(Fb,0);
        
        aP = aE + aW + aN + aS + aT + aB + (Fe-Fw) + (Fn-Fs) + (Ft-Fb);
        
        d_v(i,j,k) = alpha * An / aP;   %refer to Versteeg CFD book
        end
end

i=imax; %right
for j = 2:jmax
       %k_top=UpperBoundaryFunc(i,j,dx,dy,dz,lx,ly,lz,R,Rxy,gap,p,Ee);
       %k_top=UpperBoundaryFunc_Pp(i,j,dx,dy,dz,lx,ly,lz,R,Rxy,gap,p,Ee,u,v,w,k0,n0); %Principal stress induce deforma
        k_top=k_top1(i,j);
       for k=2:k_top-1
        De=1/2*Ae/(dx/2)*viscosity_cal(i,j,k-1,u,v,w,dx,dy,dz,k0,n0);
        Dw=1/2*Aw/(dx/2)*viscosity_cal(i-1,j,k-1,u,v,w,dx,dy,dz,k0,n0);
        Dn=An/(dy/2)*viscosity_cal(i,j,k-1,u,v,w,dx,dy,dz,k0,n0);
        Ds=As/(dy/2)*viscosity_cal(i,j-1,k-1,u,v,w,dx,dy,dz,k0,n0);
        Dt=1/2*At/(dz/2)*viscosity_cal(i,j,k,u,v,w,dx,dy,dz,k0,n0);
        Db=1/2*Ab/(dz/2)*viscosity_cal(i,j,k-1,u,v,w,dx,dy,dz,k0,n0);
        
        Fe  = 0;                                     
        Fw  = .5*rho*Ae*(u(i,j,k)+u(i,j-1,k)); 
        Fn  = .5*rho*An*(v(i,j,k)+v(i,j+1,k)); 
        Fs  = .5*rho*An*(v(i,j-1,k)+v(i,j,k));
        Ft  = .5*rho*At*(w(i,j-1,k+1)+w(i,j,k+1)); 
        Fb  = .5*rho*At*(w(i,j-1,k)+w(i,j,k));
        
        aE = 0;
        aW = Dw * A(Fw,Dw) + max(Fw,0);
        aN = Dn * A(Fn,Dn) + max(-Fn,0);
        aS = Ds * A(Fs,Ds) + max(Fs,0);
        aT = Dt * A(Ft,Dt) + max(-Ft,0);
        aB = Db * A(Fb,Db) + max(Fb,0);
        
        aP = aE + aW + aN + aS + aT + aB + (Fe-Fw) + (Fn-Fs) + (Ft-Fb);
        
        d_v(i,j,k) = alpha * An / aP;   %refer to Versteeg CFD book
        end
end


%%Apply BCs
% v_star(1:imax,1,1:kmax) = -v_star(1:imax,2,1:kmax); %south wall
% v_star(1:imax,jmax+1,1:kmax) = -v_star(1:imax,jmax,1:kmax); %north wall
% v_star(1,1:jmax+1,1:kmax) = 0.0; %left wall
% v_star(imax,1:jmax+1,1:kmax) = 0.0; %right wall

v_star(1:imax,1:jmax+1, 1) = 0.0; %bottom wall
for i=1:imax %top wall 
   for j=1:jmax+1
       %k_top=UpperBoundaryFunc(i,j,dx,dy,dz,lx,ly,lz,R,Rxy,gap,p,Ee);
       %k_top=UpperBoundaryFunc_Pp(i,j,dx,dy,dz,lx,ly,lz,R,Rxy,gap,p,Ee,u,v,w,k0,n0); %Principal stress induce deforma
       k_top=k_top1(i,j);
       [Umax, Vmax, Wmax]=UVW_boudary(i,j,lx,ly,R,Ag,wH,dx,dy);
        if k_top<kmax
            v_star(i,j,k_top) = Vmax; 
        end
   end
end

return 
end