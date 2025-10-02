function HXY=hxy(N,R2,Ee,X,Y,con_a,PXY,Redge) %kkk,scale factor
% F = @(x,y)y*sin(x)+x*cos(y);
% AA = dblquad(F,pi,2*pi,0,pi);

HXY=X*0;
H00=0.000001e-3;
for i=1:N
    for j=1:N
        
        AA=0;deltax=Redge*2/N; deltay=deltax;
        for ii=1:N
                for jj=1:N
                    rr=(X(i,j)-X(ii,jj))^2+(Y(i,j)-Y(ii,jj))^2;
                    if rr~=0
%                     if rr>=(deltax*2)^2
                       AA=AA+PXY(ii,jj)*deltax*deltay/sqrt(rr)*2/pi/Ee;
                    end
                end
        end
         %AA=0;%only for test
        if X(i,j)^2+Y(i,j)^2>con_a^2
           HXY(i,j)=H00+R2-sqrt(R2^2-X(i,j)^2-Y(i,j)^2)-(R2-sqrt(R2^2-con_a^2))+AA;
        else
            HXY(i,j)=H00+AA;
        end
    end
end
            
end