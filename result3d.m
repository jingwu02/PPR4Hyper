%the figur for the superconvergence of the method
fh=figure;
scrsz = get(0,'ScreenSize');
w=400; %width of figure
h=320;
set(fh,'Position',[scrsz(3)/2-w/2,scrsz(4)/2-w/2,w,h])

% H=[0.4 0.2 0.1 0.05 0.025];
% H2=[1.178902e-01 5.894510e-02 2.947255e-02 1.473627e-02 7.368137e-03 3.684069e-03];
% T=[0.16,0.08,0.04,0.02,0.01];

% %%%---ue=sin(t)sin(x)sin(y)sin(z)---%%%
% H=[0.3971 0.1925 0.1014 0.0525 0.0281 0.0138];
% H=[0.1925 0.1014 0.0525 0.0281 0.0138];
H=[0.2493 0.1106 0.0625 0.0269 0.0138];%sqrt(h)=t
T=[1/2 1/3 1/4 1/6 1/9];

% % k=0.01
% errfem1=[5.6142e-01 2.6353e-01 1.1900e-01 5.7239e-02 ];
% errrec1=[1.2051e+00 2.7367e-01 6.8888e-02 1.7291e-02 ];
% loglog(1./H,errfem1,'.-','linewidth',2);
% hold on;
% loglog(1./H,errfem1(end).*(H./H(end)).^(1),':','linewidth',2);
% loglog(1./H,errrec1,'.-','linewidth',2);
% loglog(1./H,errrec1(end).*(H./H(end)).^(2),':','linewidth',2);
% legend('$H^1$-error for FEM','$h^1$ convergence','$L^2$-error for recovered gradient','$h^2$ convergence','Location','southwest','Interpreter','latex')
% %title('Crank-Nicolson; Triangulation with isosceles rectangular triangles; k=0.01')
% xlabel('$1/h$','Interpreter','latex')
% ylabel('errors')

% % h=0.0357
% errfem1=[1.9539e-01 7.2118e-02 3.8963e-02 3.6508e-02 3.6445e-02];
% errrec1=[1.8792e-01 5.8598e-02 1.1162e-02 5.0868e-03 6.4655e-03];
% loglog(1./T,errfem1,'.-','linewidth',2);
% hold on;
% loglog(1./T,errfem1(end).*(T./T(end)).^(1),':','linewidth',2);
% loglog(1./T,errrec1,'.-','linewidth',2);
% loglog(1./T,errrec1(end).*(T./T(end)).^(2),':','linewidth',2);
% legend('$H^1$-error for FEM','$\tau^2$ convergence','$L^2$-error for recovered gradient','$\tau^{2}$ convergence','Location','southwest','Interpreter','latex')
% %title('Crank-Nicolson; Triangulation with isosceles rectangular triangles; n=320')
% xlabel('$1/\tau$','Interpreter','latex');
% ylabel('errors');
% %axis([2 32 4e-03 2])


% %  \tau=h
% % errfem1=[9.3850e-01 5.4165e-01 2.6085e-01 1.1863e-01 5.7195e-02 2.8055e-02 1.3900e-02];
% % errrec1=[2.6002e+00 1.1929e+00 2.6288e-01 6.6051e-02 1.6705e-02 4.1792e-03 1.0570e-03];
% errfem1=[5.2935e-01 2.5803e-01 1.1824e-01 5.7141e-02 2.8048e-02 1.3900e-02];
% errrec1=[1.1838e+00 2.4492e-01 6.1472e-02 1.5556e-02 3.8579e-03 1.0570e-03];
% loglog(1./H,errfem1,'.-','linewidth',2);
% hold on;
% loglog(1./H,errfem1(end).*(H./H(end)).^(1),':','linewidth',2);
% loglog(1./H,errrec1,'.-','linewidth',2);
% loglog(1./H,errrec1(end).*(H./H(end)).^(2),':','linewidth',2);
% legend('$H^1$-error for FEM','$h^1$ convergence','$L^2$-error for recovered gradient','$h^{2}$ convergence','Location','southwest','Interpreter','latex')
% %title('Crank-Nicolson; Triangulation with isosceles rectangular triangles; k=0.01')
% xlabel('$1/h$','Interpreter','latex')
% ylabel('errors')

%  \tau=sqrt(h)
% errfem1=[9.3850e-01 5.4165e-01 2.6085e-01 1.1863e-01 5.7195e-02 2.8055e-02 1.3900e-02];
% errrec1=[2.6002e+00 1.1929e+00 2.6288e-01 6.6051e-02 1.6705e-02 4.1792e-03 1.0570e-03];
errfem1=[3.3626e-01 1.5645e-01 8.8240e-02 3.8112e-02 1.8191e-02];
errrec1=[3.3231e-01 6.0024e-02 4.5901e-02 2.4688e-02 1.1221e-02];
loglog(1./H,errfem1,'.-','linewidth',2);
hold on;
loglog(1./H,errfem1(end).*(H./H(end)).^(1),':','linewidth',2);
loglog(1./H,errrec1,'.-','linewidth',2);
loglog(1./H,errrec1(end).*(H./H(end)).^(2),':','linewidth',2);
legend('$H^1$-error for FEM','$h^1$ convergence','$L^2$-error for recovered gradient','$h^{2}$ convergence','Location','southwest','Interpreter','latex')
%title('Crank-Nicolson; Triangulation with isosceles rectangular triangles; k=0.01')
xlabel('$1/h$','Interpreter','latex')
ylabel('errors')


% % %IRT mesh,k=0.0005(ue=sin(t)sin(x)sin(y)sin(z))
% errfem1=[4.8518e-01 2.0595e-01 9.5027e-02 4.5362e-02 ];
% errrec1=[2.3277e+00 4.2289e-01 1.8349e-01 8.1948e-02 ];
% loglog(1./H,errfem1,'.-','linewidth',2);
% hold on;
% loglog(1./H,errfem1(end).*(H./H(end)).^(1),':','linewidth',2);
% loglog(1./H,errrec1,'.-','linewidth',2);
% loglog(1./H,errrec1(end).*(H./H(end)).^(2),':','linewidth',2);
% legend('$H^1$-error for FEM','$h^1$ convergence','$L^2$-error for recovered gradient','$h^2$ convergence','Location','southwest','Interpreter','latex')
% %title('Crank-Nicolson; Triangulation with isosceles rectangular triangles; k=0.01')
% xlabel('$1/h$','Interpreter','latex')
% ylabel('errors')

% %IRT mesh,h=0.025
% T=[0.008 0.004 0.002 0.001 0.0005];
% errfem1=[4.7554e-01 7.9955e-02 4.5316e-02 4.5362e-02 4.5362e-02];
% %errrec1=[1.02E-01 2.81E-02 6.96E-03 1.68E-03 3.04E-04 8.1948e-02 ];
% loglog(1./T,errfem1,'.-','linewidth',2);
% hold on;
% loglog(1./T,errfem1(end).*(T./T(end)).^(-2/3),':','linewidth',2);
% %loglog(1./T,errfem1(end).*(T./T(end)).^(-0),':','linewidth',2);
% %loglog(1./T,errrec1,'.-','linewidth',2);
% %loglog(1./T,errrec1(end).*(T./T(end)).^(-2.2/3),':','linewidth',2);
% legend('$H^1$-error for FEM','$\tau^2$ convergence','Location','southeast','Interpreter','latex');%'$L^2$-error for recovered gradient','$\tau^{2.2}$ convergence','Location','southeast','Interpreter','latex')
% %title('Crank-Nicolson; Triangulation with isosceles rectangular triangles; n=320')
% xlabel('$\tau$','Interpreter','latex');
% ylabel('errors')

% %IRT mesh,h=0.025
% T=[0.2 0.1 0.05 0.025 0.0125];
% errfem1=[4.7554e-01 7.9955e-02 4.5316e-02 4.5659e-02 4.5316e-02];
% %errrec1=[1.02E-01 2.81E-02 6.96E-03 1.68E-03 3.04E-04 7.13E-05];
% loglog(1./T,errfem1,'.-','linewidth',2);
% hold on;
% loglog(1./T,errfem1(end).*(T./T(end)).^(-2/3),':','linewidth',2);
% %loglog(1./T,errfem1(end).*(T./T(end)).^(-0),':','linewidth',2);
% %loglog(1./T,errrec1,'.-','linewidth',2);
% %loglog(1./T,errrec1(end).*(T./T(end)).^(-2.2/3),':','linewidth',2);
% legend('$H^1$-error for FEM','$\tau^2$ convergence','Location','southeast','Interpreter','latex');%'$L^2$-error for recovered gradient','$\tau^{2.2}$ convergence','Location','southeast','Interpreter','latex')
% %title('Crank-Nicolson; Triangulation with isosceles rectangular triangles; n=320')
% xlabel('$\tau$','Interpreter','latex');
% ylabel('errors')