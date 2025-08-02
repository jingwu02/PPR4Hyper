%Mesh

% j=180;
% j=1;
 j=45;
%dt=triangulation(t,p);

%PPR nodes for the j-th node
np=size(p,2);
u=zeros(size(p,2),1);
%2d
C=sparse(t([1 1 2],:),t([2 3 3],:),1,np,np);
C=C+C.';

Gu=zeros(np,2);

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

%Plot PPR nodes and local mesh
%V=vertexAttachments(dt,[j;ij1]);

%%j=160,180

fh=figure;
scrsz = get(0,'ScreenSize');
w=400; %width of figure
h=320;
set(fh,'Position',[scrsz(3)/2-w/2,scrsz(4)/2-w/2,w,h])
%tetramesh(t(it,:),p,'facealpha',0.05);

pj(:,1)=[];
hold on;
plot(pj(1,:),pj(2,:),'b.','markersize',30);
plot(p(1,j),p(2,j),'r.','markersize',30);
pdemesh(p,e,t);
axis off;
axis tight;

