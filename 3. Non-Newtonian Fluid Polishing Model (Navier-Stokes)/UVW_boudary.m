function [Umax, Vmax, Wmax]=UVW_boudary(i,j,lx,ly,R,Ag,wH,dx,dy)

    x=(i-1)*dx; 
    y=(j-1)*dy;
    
    Umax=sqrt(R^2-(x-lx/2).^2-(y-ly/2).^2)*sin(Ag)*wH-(y-ly/2)*cos(Ag)*wH;
    Vmax=(x-lx/2)*cos(Ag)*wH;
    Wmax=(x-lx/2)*sin(Ag)*wH;

end