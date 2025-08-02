model = createpde();
model.Geometry = multicuboid(1, 1, 1);

% f=@(x,y,z,t) (3*pi^2+1)*exp(t).*sin(pi*x).*sin(pi*y).*sin(pi*z);
ff = @(location,state)location.y.^2.*tanh(location.z)/1000.*state.time;
% 波动方程 ∂²u/∂t² - ∇²u = 0
specifyCoefficients(model, 'm', 1, 'd', 0, 'c', 1, 'a', 0, 'f', ff);

% Dirichlet 边界条件
applyBoundaryCondition(model, 'dirichlet', 'Face', 1:6, 'u', 0);

% 生成网格
generateMesh(model, 'Hmax', 0.2);

% 组装矩阵
tlist=0:0.1:1;
state.time=tlist(end);
FEM = assembleFEMatrices(model,state);
K = FEM.K; M = FEM.M; F = FEM.F;