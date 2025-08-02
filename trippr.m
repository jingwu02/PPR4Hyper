function Gu=trippr(p,t,uh)
% Compute the recovered gradient by PPR
%   p: Point matrix of 2 rows
%   t: Tetrahedron matrix of 3 rows
%   u: The solution column vector.
% d=2;
% for d=2

np=size(p,2);
% I=zeros(np,1);
% I(ip)=1;

% Connection matrix
C=sparse(t([1 1 2],:),t([2 3 3],:),1,np,np);
C=C+C.';

Gu=zeros(np,2);
for j=1:np
    Cj=C(:,j);
    ij=find(Cj);
    Hij=sum((p(:,ij)-p(:,j)).^2);
    [unused,m]=sort(Hij);
    ij=ij(m);
    nj=nnz(Cj);
    nj1=nj;
    detB=-1;
    k=1;
    while detB<1e-8 %&& k<=nj
        if nj1>5
            ij1=find(Cj);
            pj=p(:,[j;ij1]);
            uj=uh([j;ij1]);
            [Guj,detB]=ppri2d(pj(1,:),pj(2,:),uj);
        end
        if nj1<=5 || detB<1e-8
            Cjk=find(C(:,ij(k)));
            Cj=Cj+sparse(Cjk,1,1,np,1);
            nj1=nnz(Cj);
            k=k+1;
        end
    end
    Gu(j,:)=Guj;
end
% 
% end

% for d=3
%   np=size(p,1); 
%   u=zeros(size(p,1),1);
%   C=sparse(t(:,[1 1 1 2 2 3]),t(:,[2 3 4 3 4 4]),1,np,np);
%   C=C+C.';
%    for j=1:np
%     Cj=C(:,j);
%     ij=find(Cj);
%     Hij=sum((p(ij,:).'-p(j,:).').^2);
%     [unused,m]=sort(Hij);
%     ij=ij(m);
%     nj=nnz(Cj);
%     nj1=nj;
%     detB=-1;
%     k=1;
%     while detB<1e-8 && k<=nj
%         if nj1>8
%             ij1=find(Cj);
%             pj=p([j;ij1],:);
%             uj=u([j;ij1]);
%             [Guj,detB]=ppri3d(pj(:,1),pj(:,2),pj(:,3),uj);
%         end
%         if nj1<=8 || detB<1e-8
%             Cjk=find(C(:,ij(k)));
%             Cj=Cj+sparse(Cjk,1,1,np,1);
%             nj1=nnz(Cj);
%             k=k+1;
%         end
%     end
%    end
% end
