function [pressure,Pp] = pres_correct(imax,jmax,kmax,dx,dy,dz,lx,ly,lz,R,Rxy,gap,rhsp,rho,d_u,d_v,d_w,p,alpha)
pressure = p;                                                               %   p = Pressure
Pp = zeros(imax,jmax,kmax);                                                 %   pressure correction 
Pp_last=Pp;

Ae=dy*dz; An=dx*dz; At=dx*dy;

while(1)
for i=1:imax
    for j=1:jmax
        k_top=UpperBoundaryFunc(i,j,dx,dy,dz,lx,ly,lz,R,Rxy,gap);
        for k=1:k_top
            if i>=2 && (i<=imax-1) && (j>=2) && (j<=jmax-1) && (k>=2) && (k<=k_top-1)              
               be=rho*d_u(i+1,j,k)*Ae; bw=rho*d_u(i,j,k)*Ae;
               cn=rho*d_v(i,j+1,k)*An; cs=rho*d_v(i,j,k)*An;
               dt=rho*d_w(i,j,k+1)*At; db=rho*d_w(i,j,k)*At; 
               a=be+bw+cn+cs+dt+db;
               Pp(i,j,k)=(be*Pp(i+1,j,k)+bw*Pp(i-1,j,k)+cn*Pp(i,j+1,k)+cs*Pp(i,j-1,k)+dt*Pp(i,j,k+1)+db*Pp(i,j,k-1)+rhsp(i,j,k))/a;
               continue;
            end

            if i==1 && j==1 && k==1             
               bw=0; be=rho*d_u(i+1,j,k)*Ae; 
               cs=0; cn=rho*d_v(i,j+1,k)*An; 
               db=0; dt=rho*d_w(i,j,k+1)*At; 
               a=be+bw+cn+cs+dt+db;
               pw=0; pe=Pp(i+1,j,k);
               ps=0; pn=Pp(i,j+1,k); 
               pb=0; pt=Pp(i,j,k+1); 
               Pp(i,j,k)=(be*pe+bw*pw+cn*pn+cs*ps+dt*pt+db*pb+rhsp(i,j,k))/a;
               continue;
            end

            if i==1 && j==1 && k==k_top             
               bw=0; be=rho*d_u(i+1,j,k)*Ae; 
               cs=0; cn=rho*d_v(i,j+1,k)*An; 
               db=rho*d_w(i,j,k-1)*At; dt=0; 
               a=be+bw+cn+cs+dt+db;
               pw=0; pe=Pp(i+1,j,k);
               ps=0; pn=Pp(i,j+1,k); 
               pb=Pp(i,j,k-1); pt=0; 
               Pp(i,j,k)=(be*pe+bw*pw+cn*pn+cs*ps+dt*pt+db*pb+rhsp(i,j,k))/a;
               continue;
            end
            
            if i==imax && j==1 && k==1             
               bw=rho*d_u(i-1,j,k)*Ae; be=0; 
               cs=0; cn=rho*d_v(i,j+1,k)*An; 
               db=0; dt=rho*d_w(i,j,k+1)*At; 
               a=be+bw+cn+cs+dt+db;
               pw=Pp(i-1,j,k); pe=0;
               ps=0; pn=Pp(i,j+1,k); 
               pb=0; pt=Pp(i,j,k+1); 
               Pp(i,j,k)=(be*pe+bw*pw+cn*pn+cs*ps+dt*pt+db*pb+rhsp(i,j,k))/a;
               continue;
            end

            if i==imax && j==1 && k==k_top             
               bw=rho*d_u(i-1,j,k)*Ae; be=0; 
               cs=0; cn=rho*d_v(i,j+1,k)*An; 
               db=rho*d_w(i,j,k-1)*At; dt=0; 
               a=be+bw+cn+cs+dt+db;
               pw=Pp(i-1,j,k); pe=0;
               ps=0; pn=Pp(i,j+1,k); 
               pb=Pp(i,j,k-1); pt=0; 
               Pp(i,j,k)=(be*pe+bw*pw+cn*pn+cs*ps+dt*pt+db*pb+rhsp(i,j,k))/a;
               continue;
            end

            if i==1 && j==jmax && k==1             
               bw=0; be=rho*d_u(i+1,j,k)*Ae; 
               cs=rho*d_v(i,j-1,k)*An; cn=0; 
               db=0; dt=rho*d_w(i,j,k+1)*At; 
               a=be+bw+cn+cs+dt+db;
               pw=0; pe=Pp(i+1,j,k);
               ps=Pp(i,j-1,k); pn=0; 
               pb=0; pt=Pp(i,j,k+1); 
               Pp(i,j,k)=(be*pe+bw*pw+cn*pn+cs*ps+dt*pt+db*pb+rhsp(i,j,k))/a;
               continue;
            end

            if i==1 && j==jmax && k==k_top             
               bw=0; be=rho*d_u(i+1,j,k)*Ae; 
               cs=rho*d_v(i,j-1,k)*An; cn=0; 
               db=rho*d_w(i,j,k-1)*At; dt=0; 
               a=be+bw+cn+cs+dt+db;
               pw=0; pe=Pp(i+1,j,k);
               ps=Pp(i,j-1,k); pn=0; 
               pb=Pp(i,j,k-1); pt=0; 
               Pp(i,j,k)=(be*pe+bw*pw+cn*pn+cs*ps+dt*pt+db*pb+rhsp(i,j,k))/a;
               continue;
            end

            if i==imax && j==jmax && k==1             
               bw=rho*d_u(i-1,j,k)*Ae; be=0; 
               cs=rho*d_v(i,j-1,k)*An; cn=0; 
               db=0; dt=rho*d_w(i,j,k+1)*At; 
               a=be+bw+cn+cs+dt+db;
               pw=Pp(i-1,j,k); pe=0;
               ps=Pp(i,j-1,k); pn=0; 
               pb=0; pt=Pp(i,j,k+1); 
               Pp(i,j,k)=(be*pe+bw*pw+cn*pn+cs*ps+dt*pt+db*pb+rhsp(i,j,k))/a;
               continue;
            end

            if i==imax && j==jmax && k==k_top             
               bw=rho*d_u(i-1,j,k)*Ae; be=0; 
               cs=rho*d_v(i,j-1,k)*An; cn=0; 
               db=rho*d_w(i,j,k-1)*At; dt=0; 
               a=be+bw+cn+cs+dt+db;
               pw=Pp(i-1,j,k); pe=0;
               ps=Pp(i,j-1,k); pn=0; 
               pb=Pp(i,j,k-1); pt=0; 
               Pp(i,j,k)=(be*pe+bw*pw+cn*pn+cs*ps+dt*pt+db*pb+rhsp(i,j,k))/a;
               continue;
            end

            if i==1 && (j~=1) && (j~=jmax) && (k~=1) && (k~=k_top)       
               bw=0; be=rho*d_u(i+1,j,k)*Ae; 
               cs=rho*d_v(i,j-1,k)*An; cn=rho*d_v(i,j+1,k)*An; 
               db=rho*d_w(i,j,k-1)*At; dt=rho*d_w(i,j,k+1)*At; 
               a=be+bw+cn+cs+dt+db;
               pw=0;
               ps=Pp(i,j-1,k); pn=Pp(i,j+1,k); 
               pb=Pp(i,j,k-1); pt=Pp(i,j,k+1); 
               Pp(i,j,k)=(be*pe+bw*pw+cn*pn+cs*ps+dt*pt+db*pb+rhsp(i,j,k))/a;
               continue;
            end

            if i==imax && (j~=1) && (j~=jmax) && (k~=1) && (k~=k_top)          
               bw=rho*d_u(i-1,j,k)*Ae; be=0; 
               cs=rho*d_v(i,j-1,k)*An; cn=rho*d_v(i,j+1,k)*An; 
               db=rho*d_w(i,j,k-1)*At; dt=rho*d_w(i,j,k+1)*At; 
               a=be+bw+cn+cs+dt+db;
               pw=Pp(i-1,j,k); pe=0;
               ps=Pp(i,j-1,k); pn=Pp(i,j+1,k); 
               pb=Pp(i,j,k-1); pt=Pp(i,j,k+1); 
               Pp(i,j,k)=(be*pe+bw*pw+cn*pn+cs*ps+dt*pt+db*pb+rhsp(i,j,k))/a;
               continue;
            end

            if j==1  && (i~=1) && (i~=imax) && (k~=1) && (k~=k_top)      
               bw=rho*d_u(i-1,j,k)*Ae; be=rho*d_u(i+1,j,k)*Ae; 
               cs=0; cn=rho*d_v(i,j+1,k)*An; 
               db=rho*d_w(i,j,k-1)*At; dt=rho*d_w(i,j,k+1)*At; 
               a=be+bw+cn+cs+dt+db;
               pw=Pp(i-1,j,k); pe=Pp(i+1,j,k);
               ps=0; pn=Pp(i,j+1,k); 
               pb=Pp(i,j,k-1); pt=Pp(i,j,k+1); 
               Pp(i,j,k)=(be*pe+bw*pw+cn*pn+cs*ps+dt*pt+db*pb+rhsp(i,j,k))/a;
               continue;
            end

            if j==jmax && (i~=1) && (i~=imax) && (k~=1) && (k~=k_top)            
               bw=rho*d_u(i-1,j,k)*Ae; be=rho*d_u(i+1,j,k)*Ae; 
               cs=rho*d_v(i,j-1,k)*An; cn=0; 
               db=rho*d_w(i,j,k-1)*At; dt=rho*d_w(i,j,k+1)*At; 
               a=be+bw+cn+cs+dt+db;
               pw=Pp(i-1,j,k); pe=Pp(i+1,j,k);
               ps=Pp(i,j-1,k); pn=0; 
               pb=Pp(i,j,k-1); pt=Pp(i,j,k+1); 
               Pp(i,j,k)=(be*pe+bw*pw+cn*pn+cs*ps+dt*pt+db*pb+rhsp(i,j,k))/a;
               continue;
            end

            if k==1 && (i~=1) && (i~=imax) && (j~=1) && (j~=jmax)             
               bw=rho*d_u(i-1,j,k)*Ae; be=rho*d_u(i+1,j,k)*Ae; 
               cs=rho*d_v(i,j-1,k)*An; cn=rho*d_v(i,j+1,k)*An; 
               db=0; dt=rho*d_w(i,j,k+1)*At; 
               a=be+bw+cn+cs+dt+db;
               pw=Pp(i-1,j,k); pe=Pp(i+1,j,k);
               ps=Pp(i,j-1,k); pn=Pp(i,j+1,k); 
               pb=0; pt=Pp(i,j,k+1); 
               Pp(i,j,k)=(be*pe+bw*pw+cn*pn+cs*ps+dt*pt+db*pb+rhsp(i,j,k))/a;
               continue;
            end

            if k==k_top && (i~=1) && (i~=imax) && (j~=1) && (j~=jmax)                
               bw=rho*d_u(i-1,j,k)*Ae; be=rho*d_u(i+1,j,k)*Ae; 
               cs=rho*d_v(i,j-1,k)*An; cn=rho*d_v(i,j+1,k)*An; 
               db=rho*d_w(i,j,k-1)*At; dt=0; 
               a=be+bw+cn+cs+dt+db;
               pw=Pp(i-1,j,k); pe=Pp(i+1,j,k);
               ps=Pp(i,j-1,k); pn=Pp(i,j+1,k); 
               pb=Pp(i,j,k-1); pt=0; 
               Pp(i,j,k)=(be*pe+bw*pw+cn*pn+cs*ps+dt*pt+db*pb+rhsp(i,j,k))/a;
               continue;
            end

        end
    end
end
Pp(1,1,1)=0;

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


