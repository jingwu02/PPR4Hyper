function [Uxc,Uyc,Uzc,vx,vy,vz,Gxc,Gyc,Gzc,vol]=ppr3d(p,t,tlist,uh,ue,uex,uey,uez)

% Ux,Uy,Uz:   grad(ue) at vertices;
% Gx,Gy,Gz:   recovered gradient at vertices;
% Uxc,Uyc,Uzc: grad(ue) at centers;
% Gxc,Gyc,Gzc: Gh(uh) at centers;
% vx,vy,vz:   grad(uh) at centers;


%PPR at t=1

%initial
% N=1:(n+1)^2;  
v=uh(:,end);

x=p(:,1);y=p(:,2);z=p(:,3);

%PPR grad
Gu=tetrappr(p,t,v);
Gx=Gu(:,1);Gy=Gu(:,2);Gz=Gu(:,3);

% %exact grad
% Ux=uex(x,y,z,tlist(end));Uy=uey(x,y,z,tlist(end));Uz=uez(x,y,z,tlist(end));
% 

%%FE grad
it1=t(:,1);it2=t(:,2);it3=t(:,3);it4=t(:,4);
ip1=p(it1,:).';ip2=p(it2,:).';ip3=p(it3,:).';ip4=p(it4,:).';
nt=size(it1,1);

dx=zeros(6,nt);
dy=dx;
dz=dx;
dx=[ip2(1,:)-ip1(1,:); ip3(1,:)-ip1(1,:); ip3(1,:)-ip2(1,:); ip4(1,:)-ip1(1,:); ip4(1,:)-ip2(1,:); ip4(1,:)-ip3(1,:)];
dy=[ip2(2,:)-ip1(2,:); ip3(2,:)-ip1(2,:); ip3(2,:)-ip2(2,:); ip4(2,:)-ip1(2,:); ip4(2,:)-ip2(2,:); ip4(2,:)-ip3(2,:)];
dz=[ip2(3,:)-ip1(3,:); ip3(3,:)-ip1(3,:); ip3(3,:)-ip2(3,:); ip4(3,:)-ip1(3,:); ip4(3,:)-ip2(3,:); ip4(3,:)-ip3(3,:)];

vol= dx(1,:).*dy(3,:).*dz(6,:)+dx(3,:).*dy(6,:).*dz(1,:)+dx(6,:).*dy(1,:).*dz(3,:)-...
 dx(1,:).*dy(6,:).*dz(3,:)-dx(3,:).*dy(1,:).*dz(6,:)-dx(6,:).*dy(3,:).*dz(1,:);

gx=[(dy(6,:).*dz(3,:)-dy(3,:).*dz(6,:))./vol;...
    (dy(2,:).*dz(6,:)-dy(6,:).*dz(2,:))./vol;...
    (dy(5,:).*dz(1,:)-dy(1,:).*dz(5,:))./vol;...
    (dy(1,:).*dz(3,:)-dy(3,:).*dz(1,:))./vol];
gy=[(dz(6,:).*dx(3,:)-dz(3,:).*dx(6,:))./vol;...
    (dz(2,:).*dx(6,:)-dz(6,:).*dx(2,:))./vol;...
    (dz(5,:).*dx(1,:)-dz(1,:).*dx(5,:))./vol;...
    (dz(1,:).*dx(3,:)-dz(3,:).*dx(1,:))./vol];
gz=[(dx(6,:).*dy(3,:)-dx(3,:).*dy(6,:))./vol;...
    (dx(2,:).*dy(6,:)-dx(6,:).*dy(2,:))./vol;...
    (dx(5,:).*dy(1,:)-dx(1,:).*dy(5,:))./vol;...
    (dx(1,:).*dy(3,:)-dx(3,:).*dy(1,:))./vol];


% gx=[(dy(6,:).*dz(3,:)-dy(3,:).*dz(6,:))./vol;...
%     (dy(6,:).*dz(2,:)-dy(2,:).*dz(6,:))./vol;...
%     (dy(5,:).*dz(1,:)-dy(1,:).*dz(5,:))./vol;...
%     (dy(3,:).*dz(1,:)-dy(1,:).*dz(3,:))./vol];
% 
% gy=[(dz(6,:).*dx(3,:)-dz(3,:).*dx(6,:))./vol;...
%     (dz(6,:).*dx(2,:)-dz(2,:).*dx(6,:))./vol;...
%     (dz(5,:).*dx(1,:)-dz(1,:).*dx(5,:))./vol;...
%     (dz(3,:).*dx(1,:)-dz(1,:).*dx(3,:))./vol];
% 
% gz=[(dx(6,:).*dy(3,:)-dx(3,:).*dy(6,:))./vol;...
%     (dx(6,:).*dy(2,:)-dx(2,:).*dy(6,:))./vol;...
%     (dx(5,:).*dy(1,:)-dx(1,:).*dy(5,:))./vol;...
%     (dx(3,:).*dy(1,:)-dx(1,:).*dy(3,:))./vol];

%Volumes of elements
vol=1/6*vol;

v1=v(it1).';
v2=v(it2).';
v3=v(it3).';
v4=v(it4).';
%FE solution at center
uhc=(v1+v2+v3+v4)/4;

%FE grad at center
vx=v1.*gx(1,:)+v2.*gx(2,:)+v3.*gx(3,:)+v4.*gx(4,:);
vy=v1.*gy(1,:)+v2.*gy(2,:)+v3.*gy(3,:)+v4.*gy(4,:);
vz=v1.*gz(1,:)+v2.*gz(2,:)+v3.*gz(3,:)+v4.*gz(4,:);

%center nodes
xc=(x(it1)+x(it2)+x(it3)+x(it4)).'/4;
yc=(y(it1)+y(it2)+y(it3)+y(it4)).'/4;
zc=(z(it1)+z(it2)+z(it3)+z(it4)).'/4;

%exact solution at center
uec=ue(xc,yc,zc,tlist(end));

%exact grads at center
Uxc=uex(xc,yc,zc,tlist(end));
Uyc=uey(xc,yc,zc,tlist(end));
Uzc=uez(xc,yc,zc,tlist(end));



%PPR on the centrial points
Gx1=Gx(it1);
Gx2=Gx(it2);
Gx3=Gx(it3);
Gx4=Gx(it4);
Gy1=Gy(it1);
Gy2=Gy(it2);
Gy3=Gy(it3);
Gy4=Gy(it4);
Gz1=Gz(it1);
Gz2=Gz(it2);
Gz3=Gz(it3);
Gz4=Gz(it4);
Gxc=(Gx1+Gx2+Gx3+Gx4).'/4;
Gyc=(Gy1+Gy2+Gy3+Gy4).'/4;
Gzc=(Gz1+Gz2+Gz3+Gz4).'/4;

%L2 error of ||u-uh||_{L2} at centers
%L2erc=sqrt(sum((uec-uhc).^2.*vol))


