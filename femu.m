%function uh=femu(b,p,e,t,tlist,c,a,f,u0,ut0)
function uh=femu(b,p,e,t,tlist,c,a,f,f0,f1)

l=length(tlist)-1;
k=tlist(2)-tlist(1);%surpose \tau=k for each time step

[K,M0,F1,O,G,H,R]=assempde(b,p,e,t,c,a,'0');%a=1 Only for the construction of the mess matrix

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


%initial values
x=p(1,:).';y=p(2,:).';

% U0=u0(x,y);U0=N'*U0;
% Ut0=ut0(x,y);Ut0=N'*Ut0;
% F=M0*f(x,y,0);F=N'*F;
% 
% U1=U0+k*Ut0+M\(k^2/2*(F-K*U0));

F0=M0*f0(x,y);F0=N'*F0;
U0=K\F0;
%Ut0=(3*pi^2)*ut0(x,y,z);Ut0=N'*Ut0;
F1=M0*f1(x,y);F1=N'*F1;%F=M0*f(x,y,z,0);F=N'*F;
U1=K\F1;

uh(:,1)=U0;
uh(:,2)=U1;
for i=2:l
    %F=f(x,y,i*k); F=M0*F; F=N'*F; uh(:,i+1)= (M+k^2*K) \ (2*M*uh(:,i)-M*uh(:,i-1)+k^2*F);%(u_tt^i,v)+a(u^i,v)=<f^{i},v>
    %F=f(x,y,(i-1)*k); F=M0*F; F=N'*F; uh(:,i+1)=2*uh(:,i)-uh(:,i-1)+M\(F-K*uh(:,i))*k^2;%(u_tt^{i+1},v)+a(u^i,v)=<f^i,v>
    F=f(x,y,(i-1)*k); F=M0*F; F=N'*F; uh(:,i+1)= (M+k^2/2*K) \ (2*M*uh(:,i)+k^2*F)-uh(:,i-1);%(u_tt^i+1,v)+a((u^{i+1}+u^{i-1})/2,v)=<f^{i},v>
end

uh=N*uh+repmat(ud,1,l+1);
