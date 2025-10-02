function [x,y]=Transform(r,theta,Rt,Qt,AlphaA) %input mm,rad, pa, mm
x=r.*cos(theta)-(Rt-Qt).*tan(AlphaA);
y=r.*sin(theta);
end