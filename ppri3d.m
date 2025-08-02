function [Gu,detB]=ppri3d(x,y,z,u)
%
J=length(x);
if J<10, fprintf('Input at least 10 points\n');end
X=x-x(1);
Y=y-y(1);
Z=z-z(1);
h=max(sqrt(X.^2+Y.^2+Z.^2));
X=X/h;
Y=Y/h;
Z=Z/h;
A=[ones(J,1) X Y Z X.^2 Y.^2 Z.^2 X.*Y Y.*Z Z.*X];
B=A.'*A;
detB=det(B);
c=B\(A.'*u);
% fprintf('determinant=%d \n',detB);
Gu=c(2:4)/h;