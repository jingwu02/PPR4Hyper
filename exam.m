
%square(pdegplot(g))
g =[ 2     2     2     2
     0     1     1     0
     1     1     0     0
     1     1     0     0
     1     0     0     1
     0     0     0     0
     1     1     1     1];
%Dirichlet boundary
b =[ 1     1     1     1
     1     1     1     1
     1     1     1     1
     1     1     1     1
     1     1     1     1
     1     1     1     1
    48    48    48    48
    48    48    48    48
    49    49    49    49
    48    48    48    48];

%mesh(����������������
% n=512;h=1/n;[p,e,t]=poimesh(g,n);%IRT mesh
[p,e,t]=initmesh(g);
for j=1:6
    [p,e,t]=refinemesh(g,p,e,t);
end
h=maxh(p,t);
n=round(1/h);

%exact solution
ue=@(x,y,t) sin(t).*sin(pi*x).*sin(pi*y);
uex=@(x,y,t) pi*sin(t).*cos(pi*x).*sin(pi*y);
uey=@(x,y,t) pi*sin(t).*sin(pi*x).*cos(pi*y);

%coefficients
k=1/500*128;
tlist=0:k:1;
f=@(x,y,t) (2*pi^2-1)*sin(t).*sin(pi*x).*sin(pi*y);
u0=@(x,y) 0*x;
ut0=@(x,y) sin(pi*x).*sin(pi*y);
f0=@(x,y) 0*x;
f1=@(x,y) k*(2*pi^2)*sin(pi*x).*sin(pi*y);
c=1;
a=1;
x=p(1,:).';
y=p(2,:).';

%solve by linear FEM
%uh = femu(b,p,e,t,tlist,c,a,f,u0,ut0);
uh = femu(b,p,e,t,tlist,c,a,f,f0,f1);

%for j=1:length(tlist), pdesurf(p,t,ue(x,y,0.01*(j-1))-uh(:,j));pause(0.2);end
%plot u at t=1, pdesurf(p,t,u(:,end))

er=ue(x,y,tlist(2))-uh(:,2);


[Ux,Uy,gx,gy,Uxc,Uyc,vx,vy,gxc,gyc,ar]=ppr(n,p,t,tlist,uh,uex,uey);

[h1er, L2erc]=error(Ux,Uy,gx,gy,Uxc,Uyc,vx,vy,gxc,gyc,ar);

% %IRT mesh
% fprintf('  h = %d, k=%3.4f : \n ||grad(u)-grad(u_h)||_H1 = %8.2e ; ||grad(u)-Gu_h||_L2 = %8.2e \n',h,tlist(2)-tlist(1),h1er,L2erc);

%Delaunay mesh
fprintf('  h = %d, k=%3.4f : \n ||grad(u)-grad(u_h)||_H1 = %8.2e ; ||grad(u)-Gu_h||_L2 = %8.2e \n',h,tlist(2)-tlist(1),h1er,L2erc);









