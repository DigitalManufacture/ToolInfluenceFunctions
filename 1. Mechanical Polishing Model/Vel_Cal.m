function Vr=Vel_Cal(x,y,AlphaA,R2,Qt,wH) %kkk,scale factor
Vrx=(R2-Qt)*sin(AlphaA)*wH-y*cos(AlphaA)*wH;
Vry=x*cos(AlphaA)*wH;
Vr=sqrt(Vrx^2+Vry^2);
end