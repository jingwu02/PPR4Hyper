function [hmax,h1er,L2erc,uh,p,t]=exam3d(h,k)
% 1. Define the Geometry of the Unit Cube
model = createpde();
L = 1; % Length of the unit cube
model.Geometry = multicuboid(L, L, L, ZOffset=-L/2); % Create a unit cube geometry


%f=@(location,state) (3*pi^2+1)*exp(state.time).*sin(pi*(location.x+0.5)).*sin(pi*(location.y+0.5)).*sin(pi*(location.z+0.5));
f=@(location,state) (3*pi^2-1)*sin(state.time).*sin(pi*(location.x+0.5)).*sin(pi*(location.y+0.5)).*sin(pi*(location.z+0.5));

%f = @(location,state)location.y.^2.*tanh(location.z)/1000;% 波动方程 ∂²u/∂t² -∇²u = f
specifyCoefficients(model, 'm', 1, 'd', 0, 'c', 1, 'a', 0, 'f', f);
% Dirichlet 边界条件
applyBoundaryCondition(model, 'dirichlet', 'Face', 1:6, 'u', 0);
% 生成网格
%h=0.4;
generateMesh(model,'GeometricOrder','linear', 'Hmax', 0.25);
pdemesh(model,'facealpha',0.05);
% axis equal
[p,e,t]=meshToPet(model.Mesh);
hmax=maxh(p,t);
%coefficients
tlist=0:k:1;
l=length(tlist)-1;
%k=tlist(2)-tlist(1);


%exact solution
% ue=@(x,y,z,t) exp(t).*sin(pi*(x+0.5)).*sin(pi*(y+0.5)).*sin(pi*z);
% uex=@(x,y,z,t) pi*exp(t).*cos(pi*(x+0.5)).*sin(pi*(y+0.5)).*sin(pi*z);
% uey=@(x,y,z,t) pi*exp(t).*sin(pi*(x+0.5)).*cos(pi*(y+0.5)).*sin(pi*z);
% uez=@(x,y,z,t) pi*exp(t).*sin(pi*(x+0.5)).*sin(pi*(y+0.5)).*cos(pi*z);
% f=@(x,y,z,t) (3*pi^2+1)*exp(t).*sin(pi*(x+0.5)).*sin(pi*(y+0.5)).*sin(pi*z);
% u0=@(x,y,z) sin(pi*(x+0.5)).*sin(pi*(y+0.5)).*sin(pi*z);
% ut0=@(x,y,z) sin(pi*(x+0.5)).*sin(pi*(y+0.5)).*sin(pi*z);
% f0=@(x,y,z) 3*pi^2*sin(pi*(x+0.5)).*sin(pi*(y+0.5)).*sin(pi*z);
% f1=@(x,y,z) 3*pi^2*(1+k+k^2/2)*sin(pi*(x+0.5)).*sin(pi*(y+0.5)).*sin(pi*z);
ue=@(x,y,z,t) sin(t).*sin(pi*(x+0.5)).*sin(pi*(y+0.5)).*sin(pi*(z+0.5));
uex=@(x,y,z,t) pi*sin(t).*cos(pi*(x+0.5)).*sin(pi*(y+0.5)).*sin(pi*(z+0.5));
uey=@(x,y,z,t) pi*sin(t).*sin(pi*(x+0.5)).*cos(pi*(y+0.5)).*sin(pi*(z+0.5));
uez=@(x,y,z,t) pi*sin(t).*sin(pi*(x+0.5)).*sin(pi*(y+0.5)).*cos(pi*(z+0.5));
f=@(x,y,z,t) (3*pi^2-1)*sin(t).*sin(pi*(x+0.5)).*sin(pi*(y+0.5)).*sin(pi*(z+0.5));
u0=@(x,y,z) 0*x;
ut0=@(x,y,z) sin(pi*(x+0.5)).*sin(pi*(y+0.5)).*sin(pi*(z+0.5));
f0=@(x,y,z) 0*x;
f1=@(x,y,z) k*3*pi^2*sin(pi*(x+0.5)).*sin(pi*(y+0.5)).*sin(pi*(z+0.5));


%解有限元
state.time=0;
% 组装矩阵
FEM = assembleFEMatrices(model,state);
K = FEM.K; M0 = FEM.M; H = FEM.H; R = FEM.R;%F = FEM.F;
%B.C.
[N,orth]=pdenullorth(H);
if size(orth,2)==0
  ud=zeros(size(K,2),1);
else
  ud=full(orth*((H*orth)\R));
end

K=N'*K*N;
M=N'*M0*N;
%initial values
x=p(1,:).';y=p(2,:).';z=p(3,:).';

% U0=u0(x,y,z);U0=N'*U0;
% Ut0=ut0(x,y,z);Ut0=N'*Ut0;
% F=M0*f(x,y,z,0);F=N'*F;
% %U1=(M+k^2/2)\(2*M*U0+F)/2+k*Ut0;
% U1=U0+k*Ut0+M\(k^2/2*(F-K*U0));
F0=M0*f0(x,y,z);F0=N'*F0;
U0=K\F0;
F1=M0*f1(x,y,z);F1=N'*F1;
U1=K\F1;

uh(:,1)=U0;
uh(:,2)=U1;
for i=2:l
    %F=f(x,y,z,i*k); F=M0*F; F=N'*F; uh(:,i+1)= (M+k^2*K) \ (2*M*uh(:,i)-M*uh(:,i-1)+k^2*F);%(u_tt^i,v)+a(u^i,v)=<f^{i},v>
    %F=f(x,y,z,(i-1)*k); F=M0*F; F=N'*F; uh(:,i+1)=2*uh(:,i)-uh(:,i-1)+M\(F-K*uh(:,i))*k^2;%(u_tt^{i+1},v)+a(u^i,v)=<f^i,v>
    F=f(x,y,z,(i-1)*k); F=M0*F; F=N'*F; uh(:,i+1)= (M+k^2/2*K) \ (2*M*uh(:,i)+k^2*F)-uh(:,i-1);%(u_tt^i+1,v)+a((u^{i+1}+u^{i-1})/2,v)=<f^{i},v>
    %ui=ue(x,y,z,i*k);ui=N'*ui;mi=max(abs(ui-uh(:,i+1)))
end

uh=N*uh+repmat(ud,1,l+1);



%for j=1:length(tlist), pdesurf(p,t,ue(x,y,0.01*(j-1))-uh(:,j));pause(0.2);end
%plot u at t=1, pdesurf(p,t,u(:,end))

er=ue(x,y,z,tlist(end))-uh(:,end);

p=p.';t=t.';

[Uxc,Uyc,Uzc,vx,vy,vz,Gxc,Gyc,Gzc,vol]=ppr3d(p,t,tlist,uh,ue,uex,uey,uez);

%H1 error of ||u-uh||_{H1} at centers
h1er=sqrt(sum(((Uxc-vx).^2+(Uyc-vy).^2+(Uzc-vz).^2).*vol));
%L2 error of ||grad(u)-Gh(uh)||_{L2} at centers
L2erc=sqrt(sum(((Uxc-Gxc).^2+(Uyc-Gyc).^2+(Uzc-Gzc).^2).*vol));


% %IRT mesh
% fprintf('  h = %d, k=%3.4f : \n ||grad(u)-grad(u_h)||_H1 = %8.2e ; ||grad(u)-Gu_h||_L2 = %8.2e \n',h,tlist(2)-tlist(1),h1er,L2erc);

%Delaunay mesh
fprintf('  h = %3.4f, k=%3.4f : \n ||u-uh||_H1 = %8.4e ; ||grad(u)-Gu_h||_L2 = %8.4e \n',hmax,tlist(2)-tlist(1),h1er,L2erc);

