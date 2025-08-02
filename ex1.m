%Solving hyperbolic equation with homogeneous Diriclet B.C.

T=1;
k=0.02;
n=T/k;

g=[    1.0000    1.0000    1.0000    1.0000
   -1.0000    0.0000    1.0000    0.0000
    0.0000    1.0000    0.0000   -1.0000
   -0.0000   -1.0000         0    1.0000
   -1.0000         0    1.0000   -0.0000
    1.0000    1.0000    1.0000    1.0000
         0         0         0         0
         0         0         0         0
         0         0         0         0
    1.0000    1.0000    1.0000    1.0000];

b=[     1     1     1     1
     1     1     1     1
     1     1     1     1
     1     1     1     1
     1     1     1     1
     1     1     1     1
    48    48    48    48
    48    48    48    48
    49    49    49    49
    48    48    48    48];

[p,e,t]=initmesh(g);
c='1';

%f='1';
f=@(x,y,t) sin(t)*(3+x.^2+y.^2)/4;
u0=@(x,y) 0*x;
ut0=@(x,y) (1-x.^2-y.^2)/4;
[K,M0,F1,O,G,H,R]=assempde(b,p,e,t,c,'1','0');%a=1 Only for the construction of the mess matrix

%B.C.
[N,orth]=pdenullorth(H);
if size(orth,2)==0
  ud=zeros(size(K,2),1);
else
  ud=full(orth*((H*orth)\R));
end

K=N'*K*N;
M=N'*M0*N;
%F=N'*F1;
u=zeros(size(K,2),n+1);

%initial value u0
x=p(1,:);y=p(2,:);
u0=u0(x,y).';
u0=N'*u0;
ut0=ut0(x,y).';
ut0=N'*ut0;
F=f(x,y,0).';
F=M0*F;
F=N'*F;
u1=u0+k*ut0+M\(k^2/2*(F-K*u0));

u(:,1)=u0;
u(:,2)=u1;

for i=2:n
    F=f(x,y,(i-1)*k).';
    F=M0*F;
    F=N'*F;
    u(:,i+1)=2*u(:,i)-u(:,i-1)+M\(F-K*u(:,i))*k^2;
end

u=N*u+repmat(ud,1,n+1);




