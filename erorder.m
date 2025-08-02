% function s=erorder(err1,err2,N1,N2)
% s=-log(err2./err1)./log(N2./N1);

%the figur for the superconvergence of the method
fh=figure;
scrsz = get(0,'ScreenSize');
w=400; %width of figure
h=320;
set(fh,'Position',[scrsz(3)/2-w/2,scrsz(4)/2-w/2,w,h])

% H=1./[10 20 40 80 160 320];
% H2=[1.178902e-01 5.894510e-02 2.947255e-02 1.473627e-02 7.368137e-03 3.684069e-03];
% T=[0.32,0.16,0.08,0.04,0.02,0.01];
H=[1/512*64 1/512*32 1/512*16 1/512*8 1/512*4 1/512*2 1/512];
H2=[1.178902e-01 5.894510e-02 2.947255e-02 1.473627e-02 7.368137e-03 3.684069e-03 1.842034e-03];
 T=[1/500*128 1/500*64 1/500*32 1/500*16 1/500*8 1/500*4 1/160];%IRT
%T=[1/500*128 1/500*64 1/500*32 1/500*16 1/500*8 1/500*4 1/500*2];
%Crank-Nicolson
%IRT mesh,k=1/500
errfem1=[2.80e-01 1.38e-01 6.85e-02 3.42e-02 1.71e-02 8.55e-03 4.27e-03];
errrec1=[2.33e-01 6.20e-02 1.57e-02 3.95e-03 9.85e-04 2.43e-04 5.76e-05];
loglog(1./H,errfem1,'.-','linewidth',2);
hold on;
loglog(1./H,errfem1(end).*(H./H(end)).^(1),':','linewidth',2);
loglog(1./H,errrec1,'.-','linewidth',2);
loglog(1./H,errrec1(end).*(H./H(end)).^(2),':','linewidth',2);
legend('$\|u(t_N)-U^N\|_1$','$O(h)$','$\|\nabla u(t_N)- G_h U^N\|$','$O(h^2)$','Location','southwest','Interpreter','latex')
%title('Crank-Nicolson; Triangulation with isosceles rectangular triangles; k=0.01')
xlabel('$1/h$','Interpreter','latex')
ylabel('errors')

% %IRT mesh,n=512
% errfem1=[4.80e-02 1.72e-02 6.16e-03 4.41e-03 4.26e-03 4.27e-03 4.27e-03];
% errrec1=[4.79e-02 1.67e-02 4.50e-03 1.11e-03 2.35e-04 2.68e-05 2.54e-05];
% loglog(1./T,errfem1,'.-','linewidth',2);
% hold on;
% loglog(1./T,errfem1(end).*(T./T(end)).^(2),':','linewidth',2);
% % loglog(1./T,errfem1(end).*(T./T(end)).^(1),':','linewidth',2);
% loglog(1./T,errrec1,'.-','linewidth',2);
% loglog(1./T,errrec1(end).*(T./T(end)).^(2),':','linewidth',2);
% legend('$\|u(t_N)-U^N\|_1$','$O(\tau^2)$','$\|\nabla u(t_N)- G_h U^N\|$','$O(\tau^{2})$','Location','northeast','Interpreter','latex')
% %title('Crank-Nicolson; Triangulation with isosceles rectangular triangles; n=320')
% xlabel('$1/\tau$','Interpreter','latex')
% ylabel('errors')

% %Delaunay mesh,k=1/500,  N is the number of total noda points
% errfem1=[1.13e-01 5.69e-02 2.86e-02 1.43e-02 7.16e-03 3.58E-03 1.79e-03];
% errrec1=[9.82e-02 2.66e-02 6.34e-03 1.53e-03 3.72e-04 9.20e-05 2.07e-05];
% loglog(1./H2,errfem1,'.-','linewidth',2);
% hold on;
% loglog(1./H2,errfem1(end).*(H2./H2(end)).^(1),':','linewidth',2);
% loglog(1./H2,errrec1,'.-','linewidth',2);
% loglog(1./H2,errrec1(end).*(H2./H2(end)).^(2),':','linewidth',2);
% legend('$\|u(t_N)-U^N\|_1$','$O(h)$','$\|\nabla u(t_N)- G_h U^N\|$','$O(h^2)$','Location','southwest','Interpreter','latex')
% %title('Crank-Nicolson; delaunay triangulation; k=0.01')
% xlabel('$1/h$','Interpreter','latex')
% ylabel('errors')

% %Delaunay mesh, h=1.842e-03
% errfem1=[4.79e-02 1.69e-02 4.87e-03 2.13e-03 1.80e-03 1.79e-03 1.79e-03];
% errrec1=[4.79e-02 1.68e-02 4.53e-03 1.15e-03 2.72e-04 5.42e-05 1.32e-05];
% loglog(1./T,errfem1,'.-','linewidth',2);
% hold on;
% loglog(1./T,errfem1(end).*(T./T(end)).^(2),':','linewidth',2);
% loglog(1./T,errrec1,'.-','linewidth',2);
% %loglog(1./T,errrec1(end).*(T./T(end)).^(2.6),':','linewidth',2);
% legend('$\|u(t_N)-U^N\|_1$','$O(\tau^2)$','$\|\nabla u(t_N)- G_h U^N\|$','Location','northeast','Interpreter','latex')
% %title('Crank-Nicolson; delaunay triangulation; N=160385')
% xlabel('$1/\tau$','Interpreter','latex')
% ylabel('errors')

% %Backward Euler
% errfem2=[2.25E-01 1.10E-01 5.49E-02 2.74E-02 1.37E-02 6.92E-03];
% errrec2=[1.63E-01 4.37E-02 1.16E-02 3.62E-03 1.64E-03 1.15E-03];
% 
% loglog(N,errfem1,'.-','linewidth',1);
% hold on;
% loglog(N,errfem1(end).*(N(end)./N).^(1),':','linewidth',1);
% loglog(N,errrec1,'.-','linewidth',1);
% loglog(N,errrec1(end).*(N(end)./N).^(2),':','linewidth',1);
% loglog(N,errfem2,'.-','linewidth',1);
% hold on;
% loglog(N,errfem2(end).*(N(end)./N).^(1),':','linewidth',1);
% loglog(N,errrec2,'.-','linewidth',1);
% loglog(N,errrec2(end).*(N(end)./N).^(1.8),':','linewidth',1);
% title('B-E; Triangulation with isosceles rectangular triangles; k=0.01')
% xlabel('n')
% ylabel('errors')





% loglog(N,err,'-.','linewidth',1);
% loglog(N,err(end).*(N(end)./N).^(1.89/3),':','linewidth',1);
% 
% xlabel('DoF');
% % xlabel('$N$','interpreter','latex');
% ylabel('Absolute value of the Error');
% %title('$M=(0,0,1)$ on an unit cube, the analytical energy is 1/6','interpreter','latex')
