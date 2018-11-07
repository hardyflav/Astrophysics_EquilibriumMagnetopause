clear all;
clc;
p = genpath('/Users/flavien/Documents/Github/CorrectionSurface');
addpath(p)

%% Definition of the grid

ThetaMaxDeg = 120;                             % Maximum value of theta on the grid, degrees
PhiMaxDeg = 90;                                % Maximum value of phi on the grid, degrees
DeltaThetaDeg = 1/2;
DeltaPhiDeg = 1/2;
ThetaSpan = (0 : DeltaThetaDeg : (ThetaMaxDeg))*pi/180;
PhiSpan = (0 : DeltaPhiDeg: PhiMaxDeg).*pi/180;

[NbPointsGrid, NbPointsGridInside,...
NbPointsThetaInside, NbPointsPhiInside,...
ThetaSpanVector, PhiSpanVector] = GridDetails(ThetaMaxDeg, PhiMaxDeg, DeltaThetaDeg, DeltaPhiDeg);



%% Construction of the initial surface

PlotsPossibilitiesInitial = ["CurveEquator", "CurveMeridianParts", "CurveMeridian", "SurfaceInitialTot", "SurfaceInitialTotContours", "SurfaceInitialTotFieldLines", "GridAnalysis", "GridWrappedAbsolute", "GridWrappedAbsoluteFieldLines", "GridWrappedRelative", "GridWrappedRelativeFieldLines"];
% Plots = ["SurfaceInitialTotContours", "SurfaceInitialTotFieldLines", "GridWrappedAbsoluteFieldLines", "GridWrappedRelativeFieldLines"];
Plots = ["CurveEquator", "CurveMeridianParts", "CurveMeridian"];
[SurfaceInitialTot, ThetaCusp, rEquator, rMeridian] = InitialGuess(ThetaMaxDeg, PhiMaxDeg, DeltaThetaDeg, DeltaPhiDeg, Plots);




%% Correction of the initial guess - Cropping + LM + Merging
SurfaceTot = SurfaceInitialTot;
NumIterationsTop = 50;
NumIterationsBottom = 165*0+50;
[SurfaceCorrected, rCroppedCorrected, SurfaceCroppedExtrapolated] = Correction(DeltaThetaDeg, DeltaPhiDeg, ThetaCusp, ThetaMaxDeg, PhiMaxDeg, rEquator, rMeridian, SurfaceTot, NumIterationsTop, NumIterationsBottom);

PlotsPossibilitiesCorrection= ["GridCroppedAnalysis", "SurfaceCroppedCorrected", "SurfaceTotCorrected", "SurfaceTotCorrectedContours", "SurfaceTotCorrectedFieldLines", "GridWrappedAbsolute", "GridWrappedAbsoluteFieldLines", "GridWrappedRelative", "GridWrappedRelativeFieldLines"];
% PlotsCorrection = ["SurfaceTotCorrectedContours", "SurfaceTotCorrectedFieldLines", "GridWrappedAbsoluteFieldLines", "GridWrappedRelativeFieldLines"];
PlotsCorrection = ["GridWrappedRelative"];
CorrectionPlot(PlotsCorrection, DeltaThetaDeg, ThetaCusp, DeltaPhiDeg, ThetaMaxDeg, PhiMaxDeg, rEquator, rMeridian, SurfaceCorrected, rCroppedCorrected, SurfaceCroppedExtrapolated)

save('Simulations/InitialAndCorrectedSurface.mat');






%{

%% EN COURS - Fringing Field: first introduction


SurfaceFullTot = vertcat(SurfaceCorrected, SurfaceCorrected(2:end, :), SurfaceCorrected(2:end, :), SurfaceCorrected(2:end-1, :));
PhiMaxFullTot = 359;
PhiMaxDeg = PhiMaxFullTot;
SurfaceCorrected = SurfaceFullTot;

PhiSpanFull = (0 : DeltaPhiDeg: PhiMaxDeg).*pi/180;


SurfaceInterior = SurfaceCorrected(2:end-1, 2:end-1);
rRightBoundary = SurfaceCorrected(:, end);
BfPrevious = [0, 0, 0];
BIMF =  [0, 0, 0].';
rss = 1;

% Taking every C-other values on the interior grid, in both directions
C=5;
SurfaceInteriorLight = SurfaceInterior(1:C:end, 1:C:end);
rSurfaceInteriorLight = SurfaceInteriorLight(:);

Nb_points_theta = ThetaMaxDeg/DeltaThetaDeg+1;
N_theta = Nb_points_theta-2; % Number of points on theta-direction of the grid, from 1 to ThetaMaxDeg-1
Nb_points_phi= PhiMaxDeg/DeltaPhiDeg+1;    % Number of points on each direction of the grid, from 0 to ThetaMaxDeg
N_phi = Nb_points_phi-2; % Number of points on theta-direction of the grid, from 1 to PhiMaxDeg-1

if (mod(N_theta,C))>0
    NThetaLight = (N_theta-mod(N_theta,C))/C+1;
else
    NThetaLight = (N_theta-mod(N_theta,C))/C;
end

if (mod(N_phi,C))>0
    NPhiLight = (N_phi-mod(N_phi,C))/C+1;
else
    NPhiLight = (N_phi-mod(N_phi,C))/C;
end
% 
% NThetaLight = (N_theta-mod(N_theta,C))/C+1;
% NPhiLight = (N_phi-mod(N_phi,C))/C+1;

% NThetaLight = N_theta/C;
% Nb_points_phi= PhiMaxDeg/DeltaPhiDeg+1;    % Number of points on each direction of the grid, from 0 to PhiMaxDeg
% N_phi = Nb_points_phi-2; % Number of points on theta-direction of the grid, from 1 to PhiMaxDeg-1
% NPhiLight = (N_phi+1)/C;

theta_vec_1Light = repmat(DeltaThetaDeg:C*DeltaThetaDeg:ThetaMaxDeg-DeltaThetaDeg,[NPhiLight 1])*pi/180;
theta_vecLight = theta_vec_1Light(:);
phi_vec_1Light = (repmat(DeltaPhiDeg:C*DeltaPhiDeg:PhiMaxDeg-DeltaPhiDeg,[1 NThetaLight])*pi/180).';
phi_vecLight = phi_vec_1Light(:);

% For every 10 degrees, in about 2 min
% For every 5 degrees, in about 6 min
% For every 2 degrees, in about 13 min
% For every 1 degrees, >12 min...
FringingFieldGrid = cell(1, length(rSurfaceInteriorLight));
for k = 1:length(rSurfaceInteriorLight)
    Positionk = [theta_vecLight(k), phi_vecLight(k), rSurfaceInteriorLight(k)];
    FringingFieldGrid{k} = FringingField(Positionk, SurfaceInterior, BfPrevious, BIMF, C, ThetaMaxDeg, PhiMaxDeg, DeltaThetaDeg, DeltaPhiDeg, rss, rEquator, rMeridian, rRightBoundary);
end

FringingFieldGridReshape = reshape(FringingFieldGrid, [NPhiLight, NThetaLight]);


% Extracting cartesian components of the fringing field
BFringingX = cellfun( @(c)c(1), FringingFieldGridReshape );
BFringingY = cellfun( @(c)c(2), FringingFieldGridReshape );
BFringingZ = cellfun( @(c)c(3), FringingFieldGridReshape );
    

% Computing spherical components of the fringing field
BFringingr = BFringingX;
BFringingTheta = BFringingX;
BFringingPhi = BFringingX;
for phi = 1:size(BFringingX, 1)
    for theta = 1:size(BFringingX, 2)
        MatrixCartesianToSpherical = [sin(theta_vec_1Light(theta))*cos(phi_vec_1Light(phi)) cos(theta_vec_1Light(theta))*cos(phi_vec_1Light(phi)) -sin(phi_vec_1Light(phi)); sin(theta_vec_1Light(theta))*sin(phi_vec_1Light(phi)) cos(theta_vec_1Light(theta))*sin(phi_vec_1Light(phi)) cos(phi_vec_1Light(phi)); cos(theta_vec_1Light(theta)) -sin(theta_vec_1Light(theta)) 0];
        SphericalComponents = [BFringingX(phi,theta), BFringingY(phi,theta), BFringingZ(phi,theta)]*MatrixCartesianToSpherical;
        BFringingr(phi,theta) = SphericalComponents(1);
        BFringingTheta(phi,theta) = SphericalComponents(2);
        BFringingPhi(phi,theta) = SphericalComponents(3);
    end
end

PcolourScale = 0.1;
figure;
pcolor(abs(BFringingr)), shading flat;
colorbar, shading flat;
caxis([0 PcolourScale]);


% Interpolating at intermediate values: BFringingr
ThetaInterp = DeltaThetaDeg:C*DeltaThetaDeg:ThetaMaxDeg-DeltaThetaDeg;
PhiInterp = DeltaPhiDeg:C*DeltaPhiDeg:PhiMaxDeg-DeltaPhiDeg;
[X, Y] = ndgrid(PhiInterp.*pi/180, ThetaInterp.*pi/180);
ValuesInterp = BFringingr;

F = griddedInterpolant(X, Y, ValuesInterp, 'spline');
PhiSpan = (DeltaPhiDeg:DeltaPhiDeg:PhiMaxDeg-DeltaPhiDeg).*pi/180;
[XInterp, YInterp] = ndgrid(PhiSpan(1:end), ThetaSpan(2:end-1));
BFringingrInterp = F(XInterp, YInterp);

figure;
pcolor(abs(BFringingrInterp)), shading flat;
colorbar
% caxis([0 PcolourScale]);



 
% Interpolating at intermediate values: BFringingTheta
ThetaInterp = DeltaThetaDeg:C*DeltaThetaDeg:ThetaMaxDeg-DeltaThetaDeg;
PhiInterp = DeltaPhiDeg:C*DeltaPhiDeg:PhiMaxDeg-DeltaPhiDeg;
[X, Y] = ndgrid(PhiInterp.*pi/180, ThetaInterp.*pi/180);
ValuesInterp = BFringingTheta;

F = griddedInterpolant(X, Y, ValuesInterp, 'spline');
[XInterp, YInterp] = ndgrid(PhiSpan(2:end-1), ThetaSpan(2:end-1));
BFringingThetaInterp = F(XInterp, YInterp);

figure;
pcolor(abs(BFringingThetaInterp)), shading flat;
colorbar
% caxis([0 PcolourScale]);


% Interpolating at intermediate values: BFringingPhi
ThetaInterp = DeltaThetaDeg:C*DeltaThetaDeg:ThetaMaxDeg-DeltaThetaDeg;
PhiInterp = DeltaPhiDeg:C*DeltaPhiDeg:PhiMaxDeg-DeltaPhiDeg;
[X, Y] = ndgrid(PhiInterp.*pi/180, ThetaInterp.*pi/180);
ValuesInterp = BFringingPhi;

F = griddedInterpolant(X, Y, ValuesInterp, 'spline');
[XInterp, YInterp] = ndgrid(PhiSpan(2:end-1), ThetaSpan(2:end-1));
BFringingPhiInterp = F(XInterp, YInterp);

figure;
pcolor(abs(BFringingPhiInterp)), shading flat;
colorbar
% caxis([0 0.1]);















    


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% EN COURS - Fringing Field: first introduction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

SurfaceTotal = ConcatenatedSurfaceInterp;
SurfaceInterior = ConcatenatedSurfaceInterp(2:end-1, 2:end-1);

% Taking every C-other values on the interior grid, in both directions
C=2;
SurfaceInteriorLight = SurfaceInterior(1:C:end, 1:C:end);
rSurfaceInteriorLight = SurfaceInteriorLight(:);

Nb_points_theta = ThetaMaxDeg/DeltaThetaDeg+1;
N_theta = Nb_points_theta-2; % Number of points on theta-direction of the grid, from 1 to ThetaMaxDeg-1
Nb_points_phi= PhiMaxDeg/DeltaPhiDeg+1;    % Number of points on each direction of the grid, from 0 to ThetaMaxDeg
N_phi = Nb_points_phi-2; % Number of points on theta-direction of the grid, from 1 to PhiMaxDeg-1

if (mod(N_theta,C))>0
    NThetaLight = (N_theta-mod(N_theta,C))/C+1;
else
    NThetaLight = (N_theta-mod(N_theta,C))/C;
end

if (mod(N_phi,C))>0
    NPhiLight = (N_phi-mod(N_phi,C))/C+1;
else
    NPhiLight = (N_phi-mod(N_phi,C))/C;
end
% 
% NThetaLight = (N_theta-mod(N_theta,C))/C+1;
% NPhiLight = (N_phi-mod(N_phi,C))/C+1;

% NThetaLight = N_theta/C;
% Nb_points_phi= PhiMaxDeg/DeltaPhiDeg+1;    % Number of points on each direction of the grid, from 0 to PhiMaxDeg
% N_phi = Nb_points_phi-2; % Number of points on theta-direction of the grid, from 1 to PhiMaxDeg-1
% NPhiLight = (N_phi+1)/C;

theta_vec_1Light = repmat(DeltaThetaDeg:C*DeltaThetaDeg:ThetaMaxDeg-DeltaThetaDeg,[NPhiLight 1])*pi/180;
theta_vecLight = theta_vec_1Light(:);
phi_vec_1Light = (repmat(DeltaPhiDeg:C*DeltaPhiDeg:PhiMaxDeg-DeltaPhiDeg,[1 NThetaLight])*pi/180).';
phi_vecLight = phi_vec_1Light(:);

% For every 10 degrees, in about 2 min
% For every 5 degrees, in about 6 min
% For every 2 degrees, in about 13 min
% For every 1 degrees, >12 min...
FringingFieldGrid = cell(1, length(rSurfaceInteriorLight));
for k = 1:length(rSurfaceInteriorLight)
    Positionk = [theta_vecLight(k), phi_vecLight(k), rSurfaceInteriorLight(k)];
    FringingFieldGrid{k} = FringingField(Positionk, SurfaceInterior, BfPrevious, BIMF, C, ThetaMaxDeg, PhiMaxDeg, DeltaThetaDeg, DeltaPhiDeg, rss, rEquator, rMeridian, rRightBoundary);
end

FringingFieldGridReshape = reshape(FringingFieldGrid, [NPhiLight, NThetaLight]);

% Extracting cartesian components of the fringing field
BFringingX = cellfun( @(c)c(1), FringingFieldGridReshape );
BFringingY = cellfun( @(c)c(2), FringingFieldGridReshape );
BFringingZ = cellfun( @(c)c(3), FringingFieldGridReshape );
    

% Computing spherical components of the fringing field
BFringingr = BFringingX;
BFringingTheta = BFringingX;
BFringingPhi = BFringingX;
for phi = 1:size(BFringingX, 1)
    for theta = 1:size(BFringingX, 2)
        MatrixCartesianToSpherical = [sin(theta_vec_1Light(theta))*cos(phi_vec_1Light(phi)) cos(theta_vec_1Light(theta))*cos(phi_vec_1Light(phi)) -sin(phi_vec_1Light(phi)); sin(theta_vec_1Light(theta))*sin(phi_vec_1Light(phi)) cos(theta_vec_1Light(theta))*sin(phi_vec_1Light(phi)) cos(phi_vec_1Light(phi)); cos(theta_vec_1Light(theta)) -sin(theta_vec_1Light(theta)) 0];
        SphericalComponents = [BFringingX(phi,theta), BFringingY(phi,theta), BFringingZ(phi,theta)]*MatrixCartesianToSpherical;
        BFringingr(phi,theta) = SphericalComponents(1);
        BFringingTheta(phi,theta) = SphericalComponents(2);
        BFringingPhi(phi,theta) = SphericalComponents(3);
    end
end

PcolourScale = 0.2;
figure;
pcolor(abs(BFringingr)), shading flat;
colorbar, shading flat;
caxis([0 PcolourScale]);


% Interpolating at intermediate values: BFringingr
ThetaInterp = DeltaThetaDeg:C*DeltaThetaDeg:ThetaMaxDeg-DeltaThetaDeg;
PhiInterp = DeltaPhiDeg:C*DeltaPhiDeg:PhiMaxDeg-DeltaPhiDeg;
[X, Y] = ndgrid(PhiInterp.*pi/180, ThetaInterp.*pi/180);
ValuesInterp = BFringingr;

F = griddedInterpolant(X, Y, ValuesInterp, 'spline');
[XInterp, YInterp] = ndgrid(PhiSpan(2:end-1), ThetaSpan(2:end-1));
BFringingrInterp = F(XInterp, YInterp);

figure;
pcolor(abs(BFringingrInterp)), shading flat;
colorbar
caxis([0 PcolourScale]);


% Interpolating at intermediate values: BFringingTheta
ThetaInterp = DeltaThetaDeg:C*DeltaThetaDeg:ThetaMaxDeg-DeltaThetaDeg;
PhiInterp = DeltaPhiDeg:C*DeltaPhiDeg:PhiMaxDeg-DeltaPhiDeg;
[X, Y] = ndgrid(PhiInterp.*pi/180, ThetaInterp.*pi/180);
ValuesInterp = BFringingTheta;

F = griddedInterpolant(X, Y, ValuesInterp, 'spline');
[XInterp, YInterp] = ndgrid(PhiSpan(2:end-1), ThetaSpan(2:end-1));
BFringingThetaInterp = F(XInterp, YInterp);

pcolor(abs(BFringingThetaInterp)), shading flat;
colorbar
caxis([0 PcolourScale]);


% Interpolating at intermediate values: BFringingPhi
ThetaInterp = DeltaThetaDeg:C*DeltaThetaDeg:ThetaMaxDeg-DeltaThetaDeg;
PhiInterp = DeltaPhiDeg:C*DeltaPhiDeg:PhiMaxDeg-DeltaPhiDeg;
[X, Y] = ndgrid(PhiInterp.*pi/180, ThetaInterp.*pi/180);
ValuesInterp = BFringingPhi;

F = griddedInterpolant(X, Y, ValuesInterp, 'spline');
[XInterp, YInterp] = ndgrid(PhiSpan(2:end-1), ThetaSpan(2:end-1));
BFringingPhiInterp = F(XInterp, YInterp);

pcolor(abs(BFringingPhiInterp)), shading flat;
colorbar
caxis([0 0.1]);




test = [BFringingThetaInterp; BFringingPhiInterp; BFringingrInterp];
test(1)
















rSurface = SurfaceInteriorLight(:);

ThetaMaxDeg = 124;
BfPrevious = [0, 0, 0];
BIMF = [0, 0, 0].';
rRightBoundary = SurfaceTotal(:, end);
% FringingFieldGrid = zeros(1, length(rSurface));

% Taking every C-other values on the interior grid, in both directions
C=2;
DeltaThetaDegLight = DeltaThetaDeg*C;
DeltaPhiDegLight = DeltaPhiDeg*C;
SurfaceInteriorLight = SurfaceInterior(1:C:end, 1:C:end);
rSurface = SurfaceInteriorLight(:);

Nb_points_theta = ThetaMaxDeg/DeltaThetaDegLight+2;    % Number of points on each direction of the grid, from 0 to ThetaMaxDeg
N_theta = Nb_points_theta-2; % Number of points on theta-direction of the grid, from 1 to ThetaMaxDeg-1
Nb_points_phi= PhiMaxDeg/DeltaPhiDegLight+1;    % Number of points on each direction of the grid, from 0 to PhiMaxDeg
N_phi = Nb_points_phi-1; % Number of points on theta-direction of the grid, from 1 to PhiMaxDeg-1







theta_vec_1 = repmat(DeltaThetaDeg:C*DeltaThetaDegLight:Max_deg_theta,[(N_phi+1)/C 1]);
theta_vec = theta_vec_1(:)*pi/180;

phi_vec_1 = repmat(DeltaPhiDeg:C*DeltaPhiDegLight:Max_deg_phi,[1 N_theta/C]).';
phi_vec = phi_vec_1(:)*pi/180;
   

FringingFieldGrid = cell(1, length(rSurface);
for k = 1:length(rSurface)/C
    Positionk = [theta_vec(C*k), phi_vec(C*k), rSurface(C*k)];
    FringingFieldGrid{k} = FringingField(Positionk, SurfaceInterior, BfPrevious, BIMF, ThetaMaxDeg, PhiMaxDeg, DeltaThetaDeg, DeltaPhiDeg, rss, rEquator, rMeridian, rRightBoundary);
end

FringingFieldGridReshape = reshape(FringingFieldGrid, [N_phi, N_theta]);









    Nb_points_theta = ThetaMaxDeg/DeltaThetaDeg+1;    % Number of points on each direction of the grid, from 0 to ThetaMaxDeg
    N_theta = Nb_points_theta-2; % Number of points on theta-direction of the grid, from 1 to ThetaMaxDeg-1
    Nb_points_phi= PhiMaxDeg/DeltaPhiDeg+1;    % Number of points on each direction of the grid, from 0 to ThetaMaxDeg
    N_phi = Nb_points_phi-2; % Number of points on theta-direction of the grid, from 1 to ThetaMaxDeg-1

    Max_deg_theta = ThetaMaxDeg-DeltaThetaDeg;
    Max_deg_phi= PhiMaxDeg-DeltaPhiDeg;
    delta_theta = DeltaThetaDeg*(pi/180);
    delta_phi = DeltaPhiDeg*(pi/180);
    
    theta_vec_1 = repmat(DeltaThetaDeg:DeltaThetaDeg:Max_deg_theta,[N_phi 1]);
    theta_vec = theta_vec_1(:)*pi/180;

    phi_vec_1 = repmat(DeltaPhiDeg:DeltaPhiDeg:Max_deg_phi,[1 N_theta]).';
    phi_vec = phi_vec_1(:)*pi/180;

% For every 10 degrees, in about 2 min
% For every 5 degrees, in about 6 min
% For every 2 degrees, in about 13 min
% For every 1 degrees, >12 min...

%     theta_vec_1 = repmat(DeltaThetaDeg:C*DeltaThetaDeg:Max_deg_theta,[(N_phi+1)/C 1]);
%     theta_vec = theta_vec_1(:)*pi/180;
% 
%     phi_vec_1 = repmat(DeltaPhiDeg:C*DeltaPhiDeg:Max_deg_phi,[1 N_theta/C]).';
%     phi_vec = phi_vec_1(:)*pi/180;
    
FringingFieldGrid = cell(1, length(rSurface)/C);
for k = 1:length(rSurface)/C
    Positionk = [theta_vec(C*k), phi_vec(C*k), rSurface(C*k)];
    FringingFieldGrid{k} = FringingField(Positionk, SurfaceInterior, BfPrevious, BIMF, ThetaMaxDeg, PhiMaxDeg, DeltaThetaDeg, DeltaPhiDeg, rss, rEquator, rMeridian, rRightBoundary);
end

FringingFieldGridBis = reshape(FringingFieldGrid, [size(SurfaceInterior, 1)/2, 2*size(SurfaceInterior, 2)/2]);
ThetaValues = theta_vec(1:C:end);
PhiValues = phi_vec(1:C:end);

    Positionk = [theta_vec(1:C:end), phi_vec(1:C:end), rSurface(1:C:end)];


Positionk = [90*pi/180, 89*pi/180, SurfaceInterior(89, 90)];
FringingField(Positionk, SurfaceInterior, BfPrevious, BIMF, ThetaMaxDeg, PhiMaxDeg, DeltaThetaDeg, DeltaPhiDeg, rss, rEquator, rMeridian, rRightBoundary)










% 
% rGuess = SurfaceTotal(2:end-1, 2:end-1);
% % Array_grid_inside / r_guess: 
% %    - theta = delta_theta : delta_theta : Max_deg_equ - delta_theta
% %    - phi = delta_phi : delta_phi : Max_deg_phi - delta_phi
% %    - size = NbPointsGridInside = (NbPointsTheta-2) * (NbPointsPhi-2)
% 
% figure (13)
% FValuesGridInsideVector = FValuesFringingField(rGuess, rEquator, r_RightBoundary(ConcatenatedSurfaceInterpVec, NbPointsPhiInside, NbPointsThetaInside, rMeridian), rMeridian, 1, ThetaMaxDeg, PhiMaxDeg, DeltaThetaDeg, DeltaPhiDeg, NbPointsGridInside, BfPrevious, BIMF);
% FValuesGridInside = reshape(FValuesGridInsideVector,[NbPointsPhiInside, NbPointsThetaInside]);
% ErrorAbsolute = (abs(FValuesGridInside(1:1:end,1:1:end)));
% pcolor(ErrorAbsolute(1:end, 1:end)),shading flat; colorbar; 
%     axes = gca;
%     xlabel('Theta','FontSize',16,'FontWeight','bold');
%     xticks(0:10/DeltaThetaDeg:NbPointsThetaInside);
%     xticklabels({(0:10/DeltaThetaDeg:NbPointsThetaInside)*DeltaThetaDeg});
%     ylabel('Phi','FontSize',16,'FontWeight','bold');
%     yticks(0:10/DeltaPhiDeg:NbPointsPhiInside);
%     yticklabels({(0:10/DeltaPhiDeg:NbPointsPhiInside)*DeltaPhiDeg});
%     set(gca, 'FontSize', 16);
%     axes.XAxisLocation = 'origin';
%     axes.YAxisLocation = 'origin';
%     set(gcf,'color','w');
%     box on










