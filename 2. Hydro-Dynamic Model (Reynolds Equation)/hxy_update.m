function HXY=hxy_update(N,Ee,X,Y,PXY,Redge,HXY_init) %kkk,scale factor
% F = @(x,y)y*sin(x)+x*cos(y);
% AA = dblquad(F,pi,2*pi,0,pi);

HXY=X*0;

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

           HXY(i,j)=HXY_init(i,j)+AA;

    end
end
            
end