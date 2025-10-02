%% plot
disp(['Total Iterations = ',int2str(iteration)])

figure; 
kkk=[2];
for i=1:length(kkk)
subplot(1,1,i);
uu=u(2:imax,2:jmax,kkk(i));
contourf(x*1000,y*1000,uu,50, 'edgecolor','none');colormap jet
colorbar;
end

%% plot
disp(['Total Iterations = ',int2str(iteration)])

figure;
% kkk=[5 10 15 20 25 30];
kkk=[16];
for i=1:length(kkk)
% subplot(2,3,i);
subplot(1,1,i);
uup=u(2:imax,kkk(i),2:kmax);
uu=zeros(imax-1,kmax-1);
for jj=1:kmax-1
   uu(:,jj)=uup(:,:,jj);
end    
contourf(x*1000,z*1000,uu',50, 'edgecolor','none');colormap jet
colorbar;
end

%% plot vector uw plane
level=16;
uup = u(2:imax,level,2:kmax);  %select u velocity field (adjust for staggered grid)
UU=zeros(imax-1,kmax-1);
for jj=1:kmax-1
   UU(:,jj)=uup(:,:,jj);
end  
wwp = w(2:imax,level,2:kmax);  %select v velocity field (adjust for staggered grid)
WW=zeros(imax-1,kmax-1);
for jj=1:kmax-1
   WW(:,jj)=wwp(:,:,jj);
end  
[X,Z]=meshgrid(x,z);   %vector plot
figure; q=quiver(X'*1000-5,Z'*1000,UU,WW,0.4);

%% plot vector uv plane
level=2;
uup = u(2:imax,2:jmax,level);  %select u velocity field (adjust for staggered grid)
UU=uup; 
vvp = v(2:imax,2:jmax,level);  %select v velocity field (adjust for staggered grid)
VV=vvp;
[X,Y]=meshgrid(x,y);   %vector plot
figure; q=quiver(X'*1000-5,Y'*1000-5,UU,VV,1);
%% slice Pp
[X,Y,Z]=meshgrid(x,y,z);
X=X-0.01/2;Y=Y-0.01/2;
X=X*1000; Y=Y*1000; Z=Z*1000;
for ii=1:33
for jj=1:33
for kk=1:44
if Pp(ii,jj,kk)>1e5
Pp(ii,jj,kk)=NaN;
end
end
end
end
uuuu=Pp(1:32,1:32,1:43);
slice(X,Y,Z,uuuu,0.00,0.00,0.006)
%% Bottom Pp
[XB,YB]=meshgrid(x,y);
XB=XB-0.01/2;YB=YB-0.01/2;
XB=XB*1000; YB=YB*1000;
mesh(XB,YB,visc_all(1:length(XB(:,1)),1:length(XB(:,1)),2));
%% normal force
F006=sum(sum(Pp(:,:,2)))*dx*dy
F009=sum(sum(Pp(:,:,2)))*dx*dy
F015=sum(sum(Pp(:,:,2)))*dx*dy
F021=sum(sum(Pp(:,:,2)))*dx*dy
theo=[0.7 0.4731 0.27]; 
exp=[0.81 0.41 0.2];
gap=[0.06 0.09 0.15];
comb=[exp' theo'];
bar(gap,comb);
%% Removal rate vs gap
theo=[0.0221 0.0102 0.0059 0.0029]/2; %kk
exp=[0.0207 0.010976 0.0069 0.0043 ]/2;
gap=[0.06 0.09 0.12 0.15];
comb=[exp' theo'];
bar(gap,comb);
%% Central deform vs gap
theo=[17.5 13.1 10.45 8.9]; %kk
gap=[0.06 0.09 0.12 0.15];
figure;plot(gap,theo);
%% Ra vs gap
theo=[70 25 8.3 5.2 3.9]; %kk
gap=[-0.06 0.06 0.12 0.18 0.21];
figure;plot(gap,theo);
%% Removal rate vs PA
theo=[0 0.0002 0.0059  0.015]/2; 
exp=[0  0.0005 0.0069 0.016]/2;
PA=[0 10 20 30];
comb=[exp' theo'];
bar(PA,comb);
%% set figure
set(gcf,'color','w');
set(0,'DefaultTextFontname','Times New Roman')
set(0,'DefaultAxesFontname','Times New Roman')
set(0,'DefaultAxesFontsize',18);
set(0,'DefaultTextFontsize',18);
xlabel('Precess angle(deg)'); ylabel('Removal rate(mm^3/min)');  zlabel('Z(mm)'); 