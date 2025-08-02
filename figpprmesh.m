%Mesh
numberOfPDEs=1;
pdem=createpde(numberOfPDEs);
gm=importGeometry(pdem,'Tetrahedron.stl');
msh=generateMesh(pdem,'hmax',9,'GeometricOrder','linear');
p=msh.Nodes.';
t=msh.Elements.';
% j=180;
 j=160;
% j=1;
dt=triangulation(t,p);

%PPR nodes for the j-th node
np=size(p,1);
u=zeros(size(p,1),1);
C=sparse(t(:,[1 1 1 2 2 3]),t(:,[2 3 4 3 4 4]),1,np,np);
C=C+C.';
    Cj=C(:,j);
    ij=find(Cj);
    Hij=sum((p(ij,:).'-p(j,:).').^2);
    [unused,m]=sort(Hij);
    ij=ij(m);
    nj=nnz(Cj);
    nj1=nj;
    detB=-1;
    k=1;
    while detB<1e-8 && k<=nj
        if nj1>8
            ij1=find(Cj);
            pj=p([j;ij1],:);
            uj=u([j;ij1]);
            [Guj,detB]=ppr(pj(:,1),pj(:,2),uj);
            %[Guj,detB]=ppri3d(pj(:,1),pj(:,2),pj(:,3),uj);
        end
        if nj1<=8 || detB<1e-8
            Cjk=find(C(:,ij(k)));
            Cj=Cj+sparse(Cjk,1,1,np,1);
            nj1=nnz(Cj);
            k=k+1;
        end
    end

%Plot PPR nodes and local mesh
V=vertexAttachments(dt,[j;ij1]);
it=zeros(1,size(t,1));
for i=1:nj1, it(V{i})=it(V{i})+1;end;
it=find(it>2);
%%j=160,180
it=V{1};
fh=figure;
scrsz = get(0,'ScreenSize');
w=400; %width of figure
h=320;
set(fh,'Position',[scrsz(3)/2-w/2,scrsz(4)/2-w/2,w,h])
tetramesh(t(it,:),p,'facealpha',0.05);
hold on;
pj(1,:)=[];
plot3(pj(:,1),pj(:,2),pj(:,3),'b.','markersize',30);
plot3(p(j,1),p(j,2),p(j,3),'r.','markersize',30);
axis off;
axis tight;


% it=[12 122 123 124 125 126 127];
% fh=figure;
% scrsz = get(0,'ScreenSize');
% w=400; %width of figure
% h=320;
% set(fh,'Position',[scrsz(3)/2-w/2,scrsz(4)/2-w/2,w,h])
% tetramesh(t(it,:),p,'facealpha',0.05);
% hold on;
% pj(1,:)=[];
% plot3(p(j,1),p(j,2),p(j,3),'r.','markersize',30);
% a=[3 4 6];
% b=[2 5 7 8 9];
% plot3(pj(a,1),pj(a,2),pj(a,3),'b.','markersize',30);
% plot3(pj(b,1),pj(b,2),pj(b,3),'k.','markersize',30);
% axis off;
% axis tight;