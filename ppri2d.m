function [Gu,detB]=ppri2d(x,y,u)
%
J=length(x);
if J<6, fprintf('Input at least 6 points\n');end
X=x.'-x(1);
Y=y.'-y(1);

h=max(sqrt(X.^2+Y.^2));

X=X/h;
Y=Y/h;

A=[ones(J,1) X Y X.^2 Y.^2 X.*Y];
B=A.'*A;
detB=det(B);
c=B\(A.'*u);
% fprintf('determinant=%d \n',detB);
Gu=c(2:3)/h;