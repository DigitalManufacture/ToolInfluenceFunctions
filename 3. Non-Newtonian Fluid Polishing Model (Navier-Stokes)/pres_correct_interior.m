function [pressure,Pp] = pres_correct_interior(imax,jmax,kmax,dx,dy,dz,lx,ly,lz,R,Rxy,gap,rhsp,rho,d_u,d_v,d_w,p,Ee,alpha,u,v,w,k0,n0,k_top1)
pressure = p;                                                               %   p = Pressure
Pp = zeros(imax,jmax,kmax);                                                 %   pressure correction 
Pp_last=Pp;

Ae=dy*dz; An=dx*dz; At=dx*dy;

while(1)
for i=2:imax-1
    for j=2:jmax-1
        %k_top=UpperBoundaryFunc(i,j,dx,dy,dz,lx,ly,lz,R,Rxy,gap,p,Ee);
        %k_top=UpperBoundaryFunc_Pp(i,j,dx,dy,dz,lx,ly,lz,R,Rxy,gap,p,Ee,u,v,w,k0,n0); %Principal stress induce deforma
        k_top=k_top1(i,j);
        for k=2:k_top-1
            b1=rho*d_u(i+1,j,k)*Ae; b2=rho*d_u(i,j,k)*Ae;
            c1=rho*d_v(i,j+1,k)*An; c2=rho*d_v(i,j,k)*An;
            d1=rho*d_w(i,j,k+1)*At; d2=rho*d_w(i,j,k)*At; 
            a=b1+b2+c1+c2+d1+d2;
            Pp(i,j,k)=(b1*Pp(i+1,j,k)+b2*Pp(i-1,j,k)+c1*Pp(i,j+1,k)+c2*Pp(i,j-1,k)+d1*Pp(i,j,k+1)+d2*Pp(i,j,k-1)+rhsp(i,j,k))/a;
        end
    end
end

max(max(max(abs(Pp-Pp_last)./abs(Pp))));  %for test
if max(max(max(abs(Pp-Pp_last)./abs(Pp))))<0.01 %correction tolerance here
    break;
end    
Pp_last=Pp;
end

pressure = p + alpha*Pp;
       
pressure(1,1,1) = 0;

return
end


