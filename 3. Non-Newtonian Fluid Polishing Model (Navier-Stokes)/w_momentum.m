function [w_star,d_w] = w_momentum(imax,jmax,kmax,lx,ly,lz,dx,dy,dz,rho,u,v,w,Ag,R,Rxy,gap,wH,k0,n0,p,Ee,alpha,k_top1)

w_star=zeros(imax,jmax,kmax+1);
d_w=zeros(imax,jmax,kmax+1);

Ae=dy*dz; An=dx*dz; At=dx*dy;
Aw=Ae; As=An; Ab=At;

A = @(F,D)( max(0, (1-0.1 * abs(F/D))^5 ) );

%%compute w_star
for i = 2:imax-1
    for j = 2:jmax-1
       %k_top=UpperBoundaryFunc(i,j,dx,dy,dz,lx,ly,lz,R,Rxy,gap,p,Ee);
       %k_top=UpperBoundaryFunc_Pp(i,j,dx,dy,dz,lx,ly,lz,R,Rxy,gap,p,Ee,u,v,w,k0,n0); %Principal stress induce deforma
       k_top=k_top1(i,j);
       for k=2:k_top-1
        De=1/2*Ae/(dx/2)*viscosity_cal(i,j,k-1,u,v,w,dx,dy,dz,k0,n0);
        Dw=1/2*Aw/(dx/2)*viscosity_cal(i-1,j,k-1,u,v,w,dx,dy,dz,k0,n0);
        Dn=1/2*An/(dy/2)*viscosity_cal(i,j,k-1,u,v,w,dx,dy,dz,k0,n0);
        Ds=1/2*As/(dy/2)*viscosity_cal(i,j-1,k-1,u,v,w,dx,dy,dz,k0,n0);
        Dt=At/(dz/2)*viscosity_cal(i,j,k,u,v,w,dx,dy,dz,k0,n0);
        Db=Ab/(dz/2)*viscosity_cal(i,j,k-1,u,v,w,dx,dy,dz,k0,n0);
        
        Fe  = .5*rho*Ae*(u(i+1,j,k)+u(i+1,j,k));                                     
        Fw  = .5*rho*Ae*(u(i-1,j,k)+u(i,j,k)); 
        Fn  = .5*rho*An*(v(i,j+1,k)+v(i,j,k)); 
        Fs  = .5*rho*An*(v(i,j-1,k)+v(i,j,k));
        Ft  = .5*rho*At*(w(i,j,k)+w(i,j,k+1)); 
        Fb  = .5*rho*At*(w(i,j,k)+w(i,j,k-1));
        
        aE = De * A(Fe,De) + max(-Fe,0);
        aW = Dw * A(Fw,Dw) + max(Fw,0);
        aN = Dn * A(Fn,Dn) + max(-Fn,0);
        aS = Ds * A(Fs,Ds) + max(Fs,0);
        aT = Dt * A(Ft,Dt) + max(-Ft,0);
        aB = Db * A(Fb,Db) + max(Fb,0);
        
        aP = aE + aW + aN + aS + aT + aB + (Fe-Fw) + (Fn-Fs) + (Ft-Fb);
        
        pressure_term = (p(i,j,k-1)-p(i,j,k)) * At;        
       
        w_star(i,j,k) = alpha/aP * ( (aE*w(i+1,j,k)+aW*w(i-1,j,k)+aN*w(i,j+1,k)+aS*w(i,j-1,k))+aT*w(i,j,k+1)+aB*w(i,j,k-1) + pressure_term ) + (1-alpha)*w(i,j,k);
        
        d_w(i,j,k) = alpha * At / aP;   %refer to Versteeg CFD book
        end
    end
end

%%set d_u for top and bottom BCs
%%they will be later used by the pressure correction equation 
%%they should not be zero, or BCs of pressure correction will get messed up
i = 1; %left
for j=2:jmax
       %k_top=UpperBoundaryFunc(i,j,dx,dy,dz,lx,ly,lz,R,Rxy,gap,p,Ee);
       %k_top=UpperBoundaryFunc_Pp(i,j,dx,dy,dz,lx,ly,lz,R,Rxy,gap,p,Ee,u,v,w,k0,n0); %Principal stress induce deforma
        k_top=k_top1(i,j);
       for k=2:k_top-1    
        De=1/2*Ae/(dx/2)*viscosity_cal(i,j,k-1,u,v,w,dx,dy,dz,k0,n0); % 
        Dw=1/2*Aw/(dx/2)*viscosity_cal(i-1,j,k-1,u,v,w,dx,dy,dz,k0,n0);
        Dn=1/2*An/(dy/2)*viscosity_cal(i,j,k-1,u,v,w,dx,dy,dz,k0,n0);
        Ds=1/2*As/(dy/2)*viscosity_cal(i,j-1,k-1,u,v,w,dx,dy,dz,k0,n0);
        Dt=At/(dz/2)*viscosity_cal(i,j,k,u,v,w,dx,dy,dz,k0,n0);
        Db=Ab/(dz/2)*viscosity_cal(i,j,k-1,u,v,w,dx,dy,dz,k0,n0);
        
        Fe  = .5*rho*Ae*(u(i+1,j,k)+u(i+1,j,k-1));                                     
        Fw  = 0; 
        Fn  = .5*rho*An*(v(i,j+1,k)+v(i,j+1,k-1)); 
        Fs  = .5*rho*An*(v(i,j,k)+v(i,j,k-1));
        Ft  = .5*rho*At*(w(i,j,k)+w(i,j,k+1)); 
        Fb  = .5*rho*At*(w(i,j,k)+w(i,j,k-1));
        
        aE = De * A(Fe,De) + max(-Fe,0);
        aW = 0;
        aN = Dn * A(Fn,Dn) + max(-Fn,0);
        aS = Ds * A(Fs,Ds) + max(Fs,0);
        aT = Dt * A(Ft,Dt) + max(-Ft,0);
        aB = Db * A(Fb,Db) + max(Fb,0);
        
        aP = aE + aW + aN + aS + aT + aB + (Fe-Fw) + (Fn-Fs) + (Ft-Fb);
        
        d_w(i,j,k) = alpha * At / aP;
       
    end
end

i = imax; %right
for j=2:jmax
       %k_top=UpperBoundaryFunc(i,j,dx,dy,dz,lx,ly,lz,R,Rxy,gap,p,Ee);
       %k_top=UpperBoundaryFunc_Pp(i,j,dx,dy,dz,lx,ly,lz,R,Rxy,gap,p,Ee,u,v,w,k0,n0); %Principal stress induce deforma
       k_top=k_top1(i,j);
       for k=2:k_top-1        
        De=1/2*Ae/(dx/2)*viscosity_cal(i,j,k-1,u,v,w,dx,dy,dz,k0,n0); % 
        Dw=1/2*Aw/(dx/2)*viscosity_cal(i-1,j,k-1,u,v,w,dx,dy,dz,k0,n0);
        Dn=An/(dy/2)*viscosity_cal(i,j,k-1,u,v,w,dx,dy,dz,k0,n0);
        Ds=As/(dy/2)*viscosity_cal(i,j-1,k-1,u,v,w,dx,dy,dz,k0,n0);
        Dt=1/2*At/(dz/2)*viscosity_cal(i,j,k,u,v,w,dx,dy,dz,k0,n0);
        Db=1/2*Ab/(dz/2)*viscosity_cal(i,j,k-1,u,v,w,dx,dy,dz,k0,n0);
        
        Fe  = 0;                                     
        Fw  = .5*rho*Ae*(u(i,j,k)+u(i,j,k-1)); 
        Fn  = .5*rho*An*(v(i,j+1,k)+v(i,j+1,k-1)); 
        Fs  = .5*rho*An*(v(i,j,k)+v(i,j,k-1));
        Ft  = .5*rho*At*(w(i,j,k)+w(i,j,k+1)); 
        Fb  = .5*rho*At*(w(i,j,k)+w(i,j,k-1));
        
        aE = 0;
        aW = Dw * A(Fw,Dw) + max(Fw,0);
        aN = Dn * A(Fn,Dn) + max(-Fn,0);
        aS = Ds * A(Fs,Ds) + max(Fs,0);
        aT = Dt * A(Ft,Dt) + max(-Ft,0);
        aB = Db * A(Fb,Db) + max(Fb,0);
        
        aP = aE + aW + aN + aS + aT + aB + (Fe-Fw) + (Fn-Fs) + (Ft-Fb);
        
        d_w(i,j,k) = alpha * At / aP;
   end
end

j = 1; %south
for i=2:imax
       %k_top=UpperBoundaryFunc(i,j,dx,dy,dz,lx,ly,lz,R,Rxy,gap,p,Ee);
       %k_top=UpperBoundaryFunc_Pp(i,j,dx,dy,dz,lx,ly,lz,R,Rxy,gap,p,Ee,u,v,w,k0,n0); %Principal stress induce deforma
       k_top=k_top1(i,j);
       for k=2:k_top-1    
        De=1/2*Ae/(dx/2)*viscosity_cal(i,j,k-1,u,v,w,dx,dy,dz,k0,n0); % 
        Dw=1/2*Aw/(dx/2)*viscosity_cal(i-1,j,k-1,u,v,w,dx,dy,dz,k0,n0);
        Dn=1/2*An/(dy/2)*viscosity_cal(i,j,k-1,u,v,w,dx,dy,dz,k0,n0);
        Ds=1/2*As/(dy/2)*viscosity_cal(i,j-1,k-1,u,v,w,dx,dy,dz,k0,n0);
        Dt=At/(dz/2)*viscosity_cal(i,j,k,u,v,w,dx,dy,dz,k0,n0);
        Db=Ab/(dz/2)*viscosity_cal(i,j,k-1,u,v,w,dx,dy,dz,k0,n0);
                
        Fe  = .5*rho*Ae*(u(i+1,j,k)+u(i+1,j,k-1));                                     
        Fw  = .5*rho*Ae*(u(i,j,k)+u(i,j,k-1)); 
        Fn  = .5*rho*An*(v(i,j+1,k)+v(i,j+1,k-1)); 
        Fs  = 0;
        Ft  = .5*rho*At*(w(i,j,k)+w(i,j,k+1)); 
        Fb  = .5*rho*At*(w(i,j,k)+w(i,j,k-1));
        
        aE = De * A(Fe,De) + max(-Fe,0);
        aW = Dw * A(Fw,Dw) + max(Fw,0);
        aN = Dn * A(Fn,Dn) + max(-Fn,0);
        aS = 0;
        aT = Dt * A(Ft,Dt) + max(-Ft,0);
        aB = Db * A(Fb,Db) + max(Fb,0);
        
        aP = aE + aW + aN + aS + aT + aB + (Fe-Fw) + (Fn-Fs) + (Ft-Fb);
        
        d_w(i,j,k) = alpha * At / aP;
       
    end
end

j = jmax; %north
for i=2:imax
       %k_top=UpperBoundaryFunc(i,j,dx,dy,dz,lx,ly,lz,R,Rxy,gap,p,Ee);
       %k_top=UpperBoundaryFunc_Pp(i,j,dx,dy,dz,lx,ly,lz,R,Rxy,gap,p,Ee,u,v,w,k0,n0); %Principal stress induce deforma
       k_top=k_top1(i,j);
       for k=2:k_top-1        
        De=1/2*Ae/(dx/2)*viscosity_cal(i,j,k-1,u,v,w,dx,dy,dz,k0,n0); % 
        Dw=1/2*Aw/(dx/2)*viscosity_cal(i-1,j,k-1,u,v,w,dx,dy,dz,k0,n0);
        Dn=1/2*An/(dy/2)*viscosity_cal(i,j,k-1,u,v,w,dx,dy,dz,k0,n0);
        Ds=1/2*As/(dy/2)*viscosity_cal(i,j-1,k-1,u,v,w,dx,dy,dz,k0,n0);
        Dt=At/(dz/2)*viscosity_cal(i,j,k,u,v,w,dx,dy,dz,k0,n0);
        Db=Ab/(dz/2)*viscosity_cal(i,j,k-1,u,v,w,dx,dy,dz,k0,n0);
        
        Fe  = .5*rho*Ae*(u(i+1,j,k)+u(i+1,j,k-1));                                     
        Fw  = .5*rho*Ae*(u(i,j,k)+u(i,j,k-1)); 
        Fn  = 0; 
        Fs  = .5*rho*An*(v(i,j,k)+v(i,j,k-1));
        Ft  = .5*rho*At*(w(i,j,k)+w(i,j,k+1)); 
        Fb  = .5*rho*At*(w(i,j,k)+w(i,j,k-1));
        
        aE = De * A(Fe,De) + max(-Fe,0);
        aW = Dw * A(Fw,Dw) + max(Fw,0);
        aN = 0;
        aS = Ds * A(Fs,Ds) + max(Fs,0);
        aT = Dt * A(Ft,Dt) + max(-Ft,0);
        aB = Db * A(Fb,Db) + max(Fb,0);
        
        aP = aE + aW + aN + aS + aT + aB + (Fe-Fw) + (Fn-Fs) + (Ft-Fb);
        
        d_w(i,j,k) = alpha * At / aP;
   end
end

%%Apply BCs
% w_star(1:imax,1,1:kmax) = 0; %south wall
% w_star(1:imax,jmax,1:kmax) = 0; %north wall
% w_star(1,1:jmax,1:kmax) = 0.0; %left wall
% w_star(imax,1:jmax,1:kmax) = 0.0; %right wall

w_star(1:imax,1:jmax, 1) = 0.0; %bottom wall
for i=1:imax %top wall 
   for j=1:jmax
       %k_top=UpperBoundaryFunc(i,j,dx,dy,dz,lx,ly,lz,R,Rxy,gap,p,Ee);
       %k_top=UpperBoundaryFunc_Pp(i,j,dx,dy,dz,lx,ly,lz,R,Rxy,gap,p,Ee,u,v,w,k0,n0); %Principal stress induce deforma
       k_top=k_top1(i,j);
       [Umax, Vmax, Wmax]=UVW_boudary(i,j,lx,ly,R,Ag,wH,dx,dy);
       if k_top<kmax
            w_star(i,j,k_top) = Wmax; 
       end
   end
end

return 
end