function Cal_k=Cal_k_coe(H,dg,miu_w,Ew) %kkk,scale factor


c=3; %relationship betwen yield stress and hardness, elastic-plastic solid is 3

dafai=(1.1*pi/c)^3*(3*(1-miu_w)^2/4)^2;
F_thresh=dafai*H^3/Ew^2*(dg/2*1e-3)^2; %in N, threshold to start yield
a=(3*(1-miu_w^2)*F_thresh*dg/2*1e-3/4/Ew)^(1/3)*1e3; %in mm, contact radius
h_elastic=a^2/(dg/2);
Cal_k=F_thresh/((H/1e6/c)*pi*(dg)*h_elastic);

end