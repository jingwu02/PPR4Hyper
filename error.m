function [h1er,L2erc]=error(Ux,Uy,gx,gy,Uxc,Uyc,vx,vy,gxc,gyc,ar)

% Ux,Uy:   grad(ue) at vertices;
% gx,gy:   PPR gradient at vertices;
% Uxc,Uyc: grad(ue) at centers;
% gxc,gyc: Gh(uh) at centers;
% vx,vy:   grad(uh) at centers;
% have to calculate: H1 half norm: [Uxc,Uyc] & [vx,vy]
%                    L2 norm     : [Uxc,Uyc] & [gxc,gyc]

%H1error of grad(u_h) at t=1
h1er=sqrt(sum(((Uxc-vx).^2+(Uyc-vy).^2).*ar));


%L2error of the ppr of grad(u_h) at t=1
L2erc=sqrt(sum(((Uxc-gxc).^2+(Uyc-gyc).^2).*ar));

% % vertices:errors between grad(ue) and Gh(u)
% erx=Ux-gx;ery=Uy-gy;
% ermax=max(sqrt(erx.^2+ery.^2));
% L2error=sqrt(sum((erx.^2+ery.^2).*ar)/3);
