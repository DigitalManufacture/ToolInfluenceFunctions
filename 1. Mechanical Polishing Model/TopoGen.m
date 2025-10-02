la=50; lb=50;
vf_min=300/60; %30 mm/min
vf_max=300/60; 
wH_min=1000;
wH_max=1000;
Qt_min=0.1;
Qt_max=0.4;
AlphaA_min=0;
AlphaA_max=30;
feedinterval=0.2;
pointspacing=0.2;
nodeNumX=la/pointspacing;
nodeNumY=lb/feedinterval;
div=1;%surface resolution coefficent
xf=linspace(0,la,nodeNumX*div);
yf=linspace(0,lb,nodeNumY*div);
[XF YF]=meshgrid(xf,yf);
% xnode=linspace(0,la,nodeNumY);
% ynode=linspace(0,lb,nodeNumY);
% [Xnode Ynode]=meshgrid(xnode,ynode);
ZF=0*XF;
VF=repmat(linspace(vf_min,vf_max,nodeNumX),nodeNumY,1);
deltaS_x=(xf(2)-xf(1))*div;
DwF=deltaS_x./VF;
deltaS_y=(yf(2)-yf(1))*div;
wHF=repmat(linspace(wH_min,wH_max,nodeNumX),nodeNumY,1); %along X
QtF=repmat((linspace(Qt_min,Qt_max,nodeNumY))',1,nodeNumX); %'along Y
AlphaAF=repmat(linspace(AlphaA_min,AlphaA_max,nodeNumX),nodeNumY,1); %along X
deviation=0.0;
     
for i=1:nodeNumY
    for j=1:nodeNumX
        [i/(nodeNumY)*100  j/(nodeNumX)*100] %progress
        [X Y Z r1]=BP_IF(QtF(i,j)+j*deviation,wHF(i,j),AlphaAF(i,j)); %if r1 is constant, putting it afterward saves time
        k2_start=round((j*deltaS_x-r1)/(deltaS_x/div)/1.2);
        k2_end=round((j*deltaS_x+r1)/(deltaS_x/div)*1.2);
        k1_start=round((i*deltaS_y-r1)/(deltaS_y/div)/1.2);
        k1_end=round((i*deltaS_y+r1)/(deltaS_y/div)*1.2);
        
%         for k1=1:nodeNumY*div
%             for k2=1:nodeNumX*div
        for k1=k1_start:k1_end
            for k2=k2_start:k2_end
                if k1>1 && (k2>1) && (k1<nodeNumY*div) && (k2<nodeNumX*div)
                if (XF(1,k2)-j*deltaS_x)^2+(YF(k1,1)-i*deltaS_y)^2<(r1)^2 %boundry threshold
                    ZF(k1,k2)=ZF(k1,k2)+interp2(X,Y,Z,XF(1,k2)-j*deltaS_x,YF(k1,1)-i*deltaS_y)*DwF(i,j);
                end
                end
            end
        end
       
    end
end
mesh(XF,YF,ZF)