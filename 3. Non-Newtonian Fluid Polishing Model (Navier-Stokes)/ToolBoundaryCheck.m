function [XT,YT,ZT] = ToolBoundaryCheck(imax,jmax,dx,dy,dz,lx,ly,lz,R,Rxy,gap,p,Ee,u,v,w,k0,n0,k_top1)

warning off
x=[]; y=[]; z=[];

for i=2:imax-1
    for j=2:jmax-1
        %k_top=UpperBoundaryFunc(i,j,dx,dy,dz,lx,ly,lz,R,Rxy,gap,p,Ee);
        %k_top=UpperBoundaryFunc_Pp(i,j,dx,dy,dz,lx,ly,lz,R,Rxy,gap,p,Ee,u,v,w,k0,n0); %Principal stress induce deforma
        k_top=k_top1(i,j);
        x=[x;dx*i]; y=[y;dy*j]; z=[z;dz*k_top];
    end
end
xt=linspace(min(x),max(x),length(x));
yt=linspace(min(y),max(y),length(y));
[XT,YT]=meshgrid(xt,yt);
ZT=griddata(x,y,z,XT,YT);
end
