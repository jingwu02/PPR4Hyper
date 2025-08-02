function [Ux, Uy, gx, gy, Uxc, Uyc, vx, vy, gxc, gyc, ar]=ppr(n,p,t,tlist,uh,uex,uey)

% Ux,Uy,Uz:   grad(ue) at vertices;
% gx,gy,gz:   PPR gradient at vertices;
% Uxc,Uyc,Uzc: grad(ue) at centers;
% gxc,gyc,gzc: Gh(uh) at centers;
% vx,vy,vz:   grad(uh) at centers;


%PPR at t=1

%initial
  N=1:(n+1)^2;
  v=uh(:,end);


x=p(1,:).';
y=p(2,:).';
Gu=trippr(p,t,v);
gx=Gu(:,1);
gy=Gu(:,2);

%exact grad
Ux=uex(x,y,tlist(end));
Uy=uey(x,y,tlist(end));
%pdesurf(p,t,Ux-gx);


%FE grad vx=v1*g1x+v2*g2x+v3*g3x;vy=v1*g1y+v2*g2y+v3*g3y
[ar,g1x,g1y,g2x,g2y,g3x,g3y]=pdetrg(p,t);
it1=t(1,:);
it2=t(2,:);
it3=t(3,:);
v1=v(it1).';
v2=v(it2).';
v3=v(it3).';
vx=v1.*g1x+v2.*g2x+v3.*g3x;
vy=v1.*g1y+v2.*g2y+v3.*g3y;
xc=(x(it1)+x(it2)+x(it3)).'/3;
yc=(y(it1)+y(it2)+y(it3)).'/3;
Uxc=uex(xc,yc,tlist(end));
Uyc=uey(xc,yc,tlist(end));
%L2er=sqrt(sum(((Uxc-vx).^2+(Uyc-vy).^2).*ar));%L2error of grad(u_h) at t=1

%PPR on the centrial points
gx1=gx(it1);
gx2=gx(it2);
gx3=gx(it3);
gy1=gy(it1);
gy2=gy(it2);
gy3=gy(it3);
gxc=(gx1+gx2+gx3).'/3;
gyc=(gy1+gy2+gy3).'/3;
%L2erc=sqrt(sum(((Uxc-gxc).^2+(Uyc-gyc).^2).*(h^2)));


%erx
% erx=Ux-gx;ery=Uy-gy;
% ermax=max(sqrt(erx.^2+ery.^2));
% L2error=sqrt(sum(((Ux-gx).^2+(Uy-gy).^2).*(h^2)));



