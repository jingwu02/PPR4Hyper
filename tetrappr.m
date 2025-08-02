% function Gu=tetrappr(p,t,u,ip)
% % Compute the recovered gradient by PPR
% %   p: Point matrix of 3 columns
% %   t: Tetrahedron matrix of 4 columns
% %   u: The solution column vector.
% %   ip: indexes of internal nodes
% 
% np=size(p,1);
% I=zeros(np,1);
% I(ip)=1;
% 
% % Connection matrix
% C=sparse(t(:,[1 1 1 2 2 3]),t(:,[2 3 4 3 4 4]),1,np,np);
% C=C+C.';
% 
% Gu=zeros(np,3);
% for j=1:np
% %     if j==17966 %5023
% %     fprintf('%d, %d ',j, sum(p(j,:).^2));
% %     end
%     Cj=C(:,j);
%     ij=find(Cj);
%     Hij=sum((p(ij,:).'-p(j,:).').^2);
%     [unused,m]=sort(Hij);
%     ij=ij(m);
%     nj=nnz(Cj);
%     nj1=nj;
%     detB=1e-10;
%     if I(j)==1
%         NeighborsofInternalPoint=1;
%     else
%         NeighborsofInternalPoint=0;
%     end
%     k=1;
%     while (nj1<=8 || detB<1e-8 || NeighborsofInternalPoint==0) && k<=nj
%         if I(ij(k))==1
%             Cjk=find(C(:,ij(k)));
%             Cj=Cj+sparse(Cjk,1,1,np,1);
%             nj1=nnz(Cj);
%             ij1=find(Cj);
%             pj=p([j;ij1],:);
%             uj=u([j;ij1]);
%             [Guj detB]=ppri3d(pj(:,1),pj(:,2),pj(:,3),uj);
%             NeighborsofInternalPoint=1;
%         end
%         k=k+1;
%     end
%     %fprintf('nj=%d\n',size(pj,1));
%     Gu(j,:)=Guj;
% end
    

% function Gu=tetrappr(p,t,u,ip)
% % Compute the recovered gradient by PPR
% %   p: Point matrix of 3 columns
% %   t: Tetrahedron matrix of 4 columns
% %   u: The solution column vector.
% %   ip: indexes of internal nodes
% 
% np=size(p,1);
% I=zeros(np,1);
% I(ip)=1;
% 
% % Connection matrix
% C=sparse(t(:,[1 1 1 2 2 3]),t(:,[2 3 4 3 4 4]),1,np,np);
% C=C+C.';
% 
% Gu=zeros(np,3);
% for j=1:np
% %     if j==17966 %5023
% %     fprintf('%d, %d ',j, sum(p(j,:).^2));
% %     end
%     Cj=C(:,j);
%     ij=find(Cj);
%     Hij=sum((p(ij,:).'-p(j,:).').^2);
%     [unused,m]=sort(Hij);
%     ij=ij(m);
%     nj=nnz(Cj);
%     nj1=nj;
%     detB=1e-10;
%     if I(j)==1
%         NeighborsofInternalPoint=1;
%     else
%         NeighborsofInternalPoint=0;
%     end
%     k=1;
%     while (nj1<=8 || detB<1e-8 || NeighborsofInternalPoint==0) && k<=nj
%         if I(ij(k))==1
%             Cjk=find(C(:,ij(k)));
%             Cj=Cj+sparse(Cjk,1,1,np,1);
%             nj1=nnz(Cj);
%             ij1=find(Cj);
%             pj=p([j;ij1],:);
%             uj=u([j;ij1]);
%             [Guj detB]=ppri3d(pj(:,1),pj(:,2),pj(:,3),uj);
%             NeighborsofInternalPoint=1;
%         end
%         k=k+1;
%     end
%     %fprintf('nj=%d\n',size(pj,1));
%     Gu(j,:)=Guj;
% end

function Gu=tetrappr(p,t,u)%,ip)
% Compute the recovered gradient by PPR
%   p: Point matrix of 3 columns
%   t: Tetrahedron matrix of 4 columns
%   u: The solution column vector.
%   ip: indexes of internal nodes

np=size(p,1);
% I=zeros(np,1);
% I(ip)=1;

% Connection matrix
C=sparse(t(:,[1 1 1 2 2 3]),t(:,[2 3 4 3 4 4]),1,np,np);
C=C+C.';

Gu=zeros(np,3);
for j=1:np
    Cj=C(:,j);
    ij=find(Cj);
    Hij=sum((p(ij,:).'-p(j,:).').^2);
    [unused,m]=sort(Hij);
    ij=ij(m);
    nj=nnz(Cj);
    nj1=nj;
    detB=-1;
    k=1;
    while detB<1e-3 && k<=nj
        if nj1>9
            ij1=find(Cj);
            pj=p([j;ij1],:);
            uj=u([j;ij1]);
            [Guj,detB]=ppri3d(pj(:,1),pj(:,2),pj(:,3),uj);
        end
        if nj1<=9 || detB<1e-3
            Cjk=find(C(:,ij(k)));
            Cj=Cj+sparse(Cjk,1,1,np,1);
            nj1=nnz(Cj);
            k=k+1;
        end
    end
    %detB
    Gu(j,:)=Guj;
end

