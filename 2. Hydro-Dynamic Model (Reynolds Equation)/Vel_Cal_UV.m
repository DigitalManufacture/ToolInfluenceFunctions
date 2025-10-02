function [U V]=Vel_Cal_UV(X,Y,AlphaA,R2,Qt,wH) %kkk,scale factor
U=(R2-Qt)*sin(AlphaA)*wH-Y.*cos(AlphaA)*wH;
V=X.*cos(AlphaA)*wH;
end