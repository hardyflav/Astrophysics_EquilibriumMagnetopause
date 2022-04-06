
clear variables
clc;
p = genpath('/Users/flavienhardy/Documents/Git/MagnetopauseTilt');
addpath(p);


% Import MP Surface
load('Phi0_5Deg.mat', 'rSol', 'ParaGrid', 'ParaSystem')
P_CartNose2KSM = ParaSystem.P_CartNose2KSM;
P_KSM2CartNose = P_CartNose2KSM^(-1);


% Build Geometry
Surface = BuildSurface(rSol, ParaGrid, ParaSystem);
figure
    surf(Surface.X_MP, Surface.Y_MP, Surface.Z_MP)
    axis equal
    
XVec = vertcat(Surface.X_MP(:), 0);
YVec = vertcat(Surface.Y_MP(:), 0);
ZVec = vertcat(Surface.Z_MP(:), 0);

Shape = alphaShape( XVec, YVec, ZVec, 'RegionThreshold', 2 );    
Shape.Alpha = 40;
plot(Shape)

[elements,nodes] = boundaryFacets(Shape);
nodes = nodes';
elements = elements';
model = createpde();
geometryFromMesh(model,nodes,elements);
pdegplot(model,'FaceLabels','on','FaceAlpha',0.5)
        

% Define Boundary Conditions: B Field on MP Surface
Interpolants.B_DipCAN = MagField_DipCAN(rSol, ParaGrid, ParaSystem); 


% Apply Boundary Conditions: MP Surface
rWhole = vertcat(rSol, flipud( rSol(1:end-1, :)) );
PhiListWhole = linspace(-pi/2, 3*pi/2, size(rWhole, 1));
[PhiGrid, ThetaGrid] = ndgrid( PhiListWhole, ParaGrid.ThetaList );

[YSurf, ZSurf, XSurf] = sph2cart(PhiGrid, pi/2-ThetaGrid, rWhole);
XSurf_KSM = P_CartNose2KSM(1,1) .* XSurf + P_CartNose2KSM(1,2) .* YSurf + P_CartNose2KSM(1,3) .* ZSurf;
YSurf_KSM = P_CartNose2KSM(2,1) .* XSurf + P_CartNose2KSM(2,2) .* YSurf + P_CartNose2KSM(2,3) .* ZSurf;
ZSurf_KSM = P_CartNose2KSM(3,1) .* XSurf + P_CartNose2KSM(3,2) .* YSurf + P_CartNose2KSM(3,3) .* ZSurf;
[nSurf_X_KSM_Val, nSurf_Y_KSM_Val, nSurf_Z_KSM_Val] = surfnorm(XSurf_KSM, YSurf_KSM, ZSurf_KSM);

Method = 'linear';
Interpolants.nP_XInterp = griddedInterpolant(PhiGrid, ThetaGrid, nSurf_X_KSM_Val, Method);
Interpolants.nP_YInterp = griddedInterpolant(PhiGrid, ThetaGrid, nSurf_Y_KSM_Val, Method);
Interpolants.nP_ZInterp = griddedInterpolant(PhiGrid, ThetaGrid, nSurf_Z_KSM_Val, Method);

Coeffunc = @(location, state) fCoeff(location, state, Interpolants, P_KSM2CartNose);
applyBoundaryCondition( model, 'neumann', 'Face', 1, 'q', 0, 'g', @(location, state) Coeffunc(location, state) );


% Apply Boundary Conditions: Closing plane
applyBoundaryCondition(model, 'neumann', 'Face', 2, 'q', 0, 'g', 0 );


% Define and Solve PDE
specifyCoefficients(model,'m',0,'d',0,'c',1,'a',0,'f',0);

hmax = 3;
generateMesh(model,'Hmax',hmax);
figure
    pdemesh(model, 'FaceAlpha',0.5); 
    axis equal

i = 1;
C_n = findNodes(model.Mesh,'Region','Face',i);

result = solvepde(model);


% Interpolate and Visualise Result
[X,Y,Z] = meshgrid(-0:10,-10:10,-10:10);
V = interpolateSolution(result,X,Y,Z);
V = reshape(V,size(X));

figure
hold on
    colormap jet
    contourslice(X,Y,Z,V,[],[0],[])
    xlabel('x')
    ylabel('y')
    zlabel('z')
    colorbar
    axis equal
    
PlotShape = plot(Shape);
set(PlotShape, 'FaceAlpha', 0.1, 'EdgeColor', 'none')

axis([-40 40 -40 40 -40 40])




[X,Y,Z] = meshgrid(-0:10,0,-10:10);
V2 = interpolateSolution(result,X,Y,Z);
V2 = squeeze( reshape( V2,size(X) ) );


