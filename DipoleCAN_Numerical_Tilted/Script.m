clear variables
clc;
p = genpath('/Users/flavien/Documents/Projects/ResearchPhysics/Github/EquilibriumMagnetopause/DipoleCAN_Numerical_Tilted');
addpath(p);

global DegToRad
global RadToDeg
global ParaSystem
DegToRad = pi/180;
RadToDeg = 180/pi;

Inclination = 1;


%% ------------------------------------------------------------------------
% Set Up: Planet and Pressure Balance Parameters

    Planet = 'Saturn';
    ParaSystem = System_ini(Planet, Inclination);
    ParaSystem.BShielding.X = 0;
    ParaSystem.BShielding.Y = 0;
    ParaSystem.BShielding.Z = 0;

%% ------------------------------------------------------------------------
% Computing CAN Disk X-Y-Z KSM Components
%{

    rList_Rp = (10 : 1 : 35) ;
    PhiList = (-90:3:90) .' .* DegToRad;
    ThetaList = (0:3:90) .' .* DegToRad;
    [GridPhi, GridTheta, GridR] = ndgrid( PhiList, ThetaList, rList_Rp );

    BCAN = B_CAN_Corrected(CAN_DiskParameters, GridR, GridTheta, GridPhi);
        
    BCAN.XInterp = griddedInterpolant(GridPhi, GridTheta, GridR, BCAN.X);
    BCAN.YInterp = griddedInterpolant(GridPhi, GridTheta, GridR, BCAN.Y);
    BCAN.ZInterp = griddedInterpolant(GridPhi, GridTheta, GridR, BCAN.Z);

    ParaSystem.BCAN = BCAN;
    [XGrid_MB, YGrid_MB, ZGrid_MB] = sph2cart(GridPhi, pi/2-GridTheta, GridR);
    
figure;
    PlotQuiver = quiver3(ZGrid_MB, XGrid_MB, YGrid_MB, BCAN.X, BCAN.Y, BCAN.Z);
    set(PlotQuiver, 'AutoScaleFactor', 6, 'LineWidth', 1, 'MaxHeadSize', 1);
    axis equal
%}
%% ------------------------------------------------------------------------
% Computing ROTATED CAN Disk X-Y-Z KSM Components FROM KSM coordinates

ParaSystem.OnOff_CAN = 1;
CAN_DiskParameters = ParaSystem.CAN_DiskParameters;

XList_KSM = (-36 : 2 : 36) .';   
YList_KSM = (-36 : 2 : 36) .';   
ZList_KSM = (-36 : 2 : 36) .';   
[Position_KSM] = [XList_KSM, YList_KSM, ZList_KSM] .';

[XGrid_KSM, YGrid_KSM, ZGrid_KSM] = ndgrid(XList_KSM, YList_KSM, ZList_KSM);

% ----  Rotating back to dipole-aligned case: Position_Aligned = [X; Y; Z]
P_AlignedToRot = ParaSystem.Tilt.RotationPhi * ParaSystem.Tilt.RotationTheta;
P_RotToAligned = P_AlignedToRot.';

XGrid_Aligned = P_RotToAligned(1,1) .* XGrid_KSM + P_RotToAligned(1,2) .* YGrid_KSM  + P_RotToAligned(1,3) .* ZGrid_KSM;
YGrid_Aligned = P_RotToAligned(2,1) .* XGrid_KSM + P_RotToAligned(2,2) .* YGrid_KSM  + P_RotToAligned(2,3) .* ZGrid_KSM;
ZGrid_Aligned = P_RotToAligned(3,1) .* XGrid_KSM + P_RotToAligned(3,2) .* YGrid_KSM  + P_RotToAligned(3,3) .* ZGrid_KSM;

ZGrid_Aligned_MB = XGrid_Aligned;
XGrid_Aligned_MB = YGrid_Aligned;
YGrid_Aligned_MB = ZGrid_Aligned;

[Phi_Aligned, Elevation_Aligned, r_Aligned] = cart2sph( XGrid_Aligned_MB, YGrid_Aligned_MB, ZGrid_Aligned_MB );
Theta_Aligned = pi/2 - Elevation_Aligned;

Phi_Aligned = wrapToPi(Phi_Aligned);

BCAN_Aligned = B_CAN_Corrected(CAN_DiskParameters, r_Aligned, Theta_Aligned, Phi_Aligned);

BCAN_Rotated.X = P_AlignedToRot(1,1) .* BCAN_Aligned.X + P_AlignedToRot(1,2) .* BCAN_Aligned.Y  + P_AlignedToRot(1,3) .* BCAN_Aligned.Z;
BCAN_Rotated.Y = P_AlignedToRot(2,1) .* BCAN_Aligned.X + P_AlignedToRot(2,2) .* BCAN_Aligned.Y  + P_AlignedToRot(2,3) .* BCAN_Aligned.Z;
BCAN_Rotated.Z = P_AlignedToRot(3,1) .* BCAN_Aligned.X + P_AlignedToRot(3,2) .* BCAN_Aligned.Y  + P_AlignedToRot(3,3) .* BCAN_Aligned.Z;


BCAN_Rotated.XInterp = griddedInterpolant(XGrid_KSM, YGrid_KSM, ZGrid_KSM, BCAN_Rotated.X, 'spline');
BCAN_Rotated.YInterp = griddedInterpolant(XGrid_KSM, YGrid_KSM, ZGrid_KSM, BCAN_Rotated.Y, 'spline');
BCAN_Rotated.ZInterp = griddedInterpolant(XGrid_KSM, YGrid_KSM, ZGrid_KSM, BCAN_Rotated.Z, 'spline');

XList = (-20:1:30);
YList = (-20:1:20);
ZList = (-20:1:20);

[X, Y, Z] = ndgrid( XList, YList, ZList );
figure;
    PlotQuiver = quiver3(X, Y, Z, BCAN_Rotated.XInterp(X, Y, Z), BCAN_Rotated.YInterp(X, Y, Z), BCAN_Rotated.ZInterp(X, Y, Z));
    set(PlotQuiver, 'AutoScaleFactor', 2, 'LineWidth', 1);
    axis equal
    
BCAN = BCAN_Rotated;
ParaSystem.BCAN = BCAN;
    
    %% -------------------------------------------------------------------------------
% 3D Plot of the System

    RotationOptions = ["DipoleTilted", "DipoleAligned"];
    Rotation = RotationOptions(1);
    Plot = Plot3D_ini(Rotation, ParaSystem);


%% -------------------------------------------------------------------------------
% Finding Position of Anchor Points: Stagnation Point and Terminator Point

% ----  Finding B.v=0 locus in NMM meridian plane + Stagnation Point
    ParaSystem = NoseDetermination(ParaSystem);

% ----  Finding Bxv=0 locus in NMM meridian plane + Terminator Point
    ParaSystem = TerminatorPoint(ParaSystem);



%% -------------------------------------------------------------------------------
% PB in NMM, Southern Branch
    
% ----  Setting up Guess and Grid
    Nose_r0 = ParaSystem.Nose.KSMU_Rp .* ParaSystem.Rp ./ ParaSystem.r0;
    ParaSystem.NoseDirection = Nose_r0 ./ norm(Nose_r0);
    Nose_Theta = acos(ParaSystem.NoseDirection(1)) ;
    
    A1Theta = acos(ParaSystem.A1.TiltedDipole.KSM_Rp(1) ./ norm(ParaSystem.A1.TiltedDipole.KSM_Rp));
    A1R = norm(ParaSystem.A1.TiltedDipole.KSM_Rp) .* ParaSystem.Rp ./ ParaSystem.r0;
    
    ParaGrid.DeltaTheta = (1/3).* DegToRad;
    ParaGrid.PhiList = -pi/2;
    ParaGrid.ThetaList = (0 : ParaGrid.DeltaTheta * RadToDeg : (A1Theta-Nose_Theta) * RadToDeg).' .* DegToRad;
%     ParaGrid.ThetaList = (0 : ParaGrid.DeltaTheta * RadToDeg : 90 ).' .* DegToRad;
%     ParaGrid.ThetaList = ((A1Theta-Nose_Theta) * RadToDeg : -ParaGrid.DeltaTheta * RadToDeg : 80).' .* DegToRad;
    [ParaGrid.PhiGrid, ParaGrid.ThetaGrid] = ndgrid(ParaGrid.PhiList, ParaGrid.ThetaList);
    
    rGrid_Func = @(kappa, theta) norm(Nose_r0)* (2./(1+cos(theta))).^(kappa).';
    
    KappaGuess = 0.6;
    kappaSol = fsolve( @(kappa) rGrid_Func(kappa, A1Theta-Nose_Theta)  - A1R, KappaGuess );
    
    rGrid = norm(Nose_r0)* (2./(1+cos(ParaGrid.ThetaGrid))).^(kappaSol) .';
    rGrid = rGrid(1:end);
    
    ParaGrid.ThetaList = ParaGrid.ThetaList(1:end);
    [ParaGrid.PhiGrid, ParaGrid.ThetaGrid] = ndgrid(ParaGrid.PhiList, ParaGrid.ThetaList);
    
    [X, Y] = pol2cart(ParaGrid.ThetaList, rGrid);
    figure;
    hold on
    plot(Y, X)
    axis equal

% ---- Optimisation Procedure    
    MaxIter = 50;
    MaxEval = 5e5;
    ParaGrid.epsilon = 1e-6;
    PBGrid_Func = @(r) PB_Grid_NMM(r, ParaGrid, ParaSystem);
    Algorithms = ["levenberg-marquardt"; "Trust-Region-Dogleg"; "Trust-Region" ];
    opts = optimoptions('fsolve', 'Display','iter', 'Maxiter', MaxIter, 'MaxFunctionEvaluations', MaxEval, 'SpecifyObjectiveGradient', false, 'CheckGradients', false, 'Algorithm', Algorithms(1), ...
                        'StepTolerance', 1e-10);

    rGuess = rGrid(:);
    rSolNMM_1 = fsolve(PBGrid_Func, rGuess, opts);
    
    rSolNMM_1 = smoothdata(rSolNMM_1, 'movmean', 4);
    [X2, Y2] = pol2cart(ParaGrid.ThetaList, rSolNMM_1);

    hold on
    plot(Y2, X2)
    axis equal

    rGuess = rSolNMM_1;
    

% Packaging Output as Interpolated Function
rNMM_South_Interp = griddedInterpolant( vertcat(ParaGrid.ThetaList), rSolNMM_1);



    %% PB in NMM, Northern Branch

    % Setting up Guess and Grid    
        ParaSystemNorth = TerminatorPointNorth(ParaSystem);

        Nose_r0 = ParaSystemNorth.Nose.KSMU_Rp .* ParaSystem.Rp ./ ParaSystem.r0;
        ParaSystemNorth.NoseDirection = Nose_r0 ./ norm(Nose_r0);
        Nose_Theta = acos(ParaSystem.NoseDirection(1)) ;
        A1Theta = acos(ParaSystemNorth.A1.TiltedDipole.KSM_Rp(1) ./ norm(ParaSystemNorth.A1.TiltedDipole.KSM_Rp));
        A1R = norm(ParaSystemNorth.A1.TiltedDipole.KSM_Rp) .* ParaSystem.Rp ./ ParaSystem.r0;

        ParaGrid.DeltaTheta = (1/2).* DegToRad;
        ParaGrid.PhiList = pi/2;
        ParaGrid.ThetaList = (0 : ParaGrid.DeltaTheta * RadToDeg : (A1Theta+Nose_Theta) * RadToDeg).' .* DegToRad;
        [ParaGrid.PhiGrid, ParaGrid.ThetaGrid] = ndgrid(ParaGrid.PhiList, ParaGrid.ThetaList);

        rGrid_Func = @(kappa, theta) norm(Nose_r0)* (2./(1+cos(theta))).^(kappa).';

        KappaGuess = 0.6;
        kappaSol = fsolve( @(kappa) rGrid_Func(kappa, A1Theta+Nose_Theta)  - A1R, KappaGuess );

        rGrid = norm(Nose_r0)* (2./(1+cos(ParaGrid.ThetaGrid))).^(kappaSol) .';
        rGrid = rGrid(1:end);

        ParaGrid.ThetaList = ParaGrid.ThetaList(1:end);
        [ParaGrid.PhiGrid, ParaGrid.ThetaGrid] = ndgrid(ParaGrid.PhiList, ParaGrid.ThetaList);

        [X, Y] = pol2cart(ParaGrid.ThetaList, rGrid);
        figure;
        hold on
        plot(Y, X)
        axis equal

    % Optimisation Procedure    
        MaxIter = 50;
        MaxEval = 5e5;
        ParaGrid.epsilon = 1e-6;
        PBGrid_Func = @(r) PB_Grid_NMM(r, ParaGrid, ParaSystem);
        Algorithms = ["levenberg-marquardt"; "Trust-Region-Dogleg"; "Trust-Region" ];
        opts = optimoptions('fsolve', 'Display','iter', 'Maxiter', MaxIter, 'MaxFunctionEvaluations', MaxEval, 'SpecifyObjectiveGradient', false, 'CheckGradients', false, 'FiniteDifferenceType', 'central', 'Algorithm', Algorithms(1));

        rGuess = rGrid(:);
        rSolNMM_1 = fsolve(PBGrid_Func, rGuess, opts);
        rGridNMM_1 = rSolNMM_1;
        
        rSolNMM_1 = smoothdata(rSolNMM_1, 'movmean', 4);
        [X2, Y2] = pol2cart(ParaGrid.ThetaList, rSolNMM_1);

        plot(Y2, X2)
        axis equal    

    % Packaging Output as Interpolated Function    
        rNMM_North_Interp = griddedInterpolant(ParaGrid.ThetaList, rSolNMM_1);

        
     %% Building Guess Surface with two NMM boundaries
        rNMM_North = rNMM_North_Interp(ParaGrid.ThetaList);
        rNMM_South = rNMM_South_Interp(ParaGrid.ThetaList);
        rMiddle = norm(Nose_r0)* (2./(1+cos(ParaGrid.ThetaList))).^(0.42) ;
        rNMM_List = [rNMM_South, rMiddle, rNMM_North];

        PhiNMM = [-pi/2; 0; pi/2];
        [PhiGrid_NMM, ThetaGrid_NMM] = ndgrid(PhiNMM, ParaGrid.ThetaList);
        rInterp = griddedInterpolant(PhiGrid_NMM, ThetaGrid_NMM, rNMM_List.', 'spline');

        ParaGrid.DeltaPhi = 3 .* DegToRad;
        ParaGrid.DeltaTheta = 3 .* DegToRad;
        ParaGrid.ThetaList = (0: ParaGrid.DeltaTheta * RadToDeg :90).' .* DegToRad;
        ParaGrid.PhiList = (-90: ParaGrid.DeltaPhi * RadToDeg : 90).' .* DegToRad;
        [ParaGrid.PhiGrid, ParaGrid.ThetaGrid] = ndgrid(ParaGrid.PhiList, ParaGrid.ThetaList);
        rGuess = rInterp(ParaGrid.PhiGrid, ParaGrid.ThetaGrid);

        rGuess(:, end) = 2 .* rGuess(:, end-1) - rGuess(:, end-2);
        % TEST
%          rGuess(:, 3) = (1/2) .* (  rGuess(:, 1) +  rGuess(:, 2) );
%          rGuess(:, end) = (1/2) .* (  rGuess(:, 1) +  rGuess(:, 2) );
%          rGuess(end-10:end, :) = smoothdata(rGuess(end-10:end, :), 2);
        % TEST
        rGuess_ini = rGuess;


    %% -------------------------------------------------------------------------------
    % Evaluating Guess Surface
        rGrid = rGuess;

        ResidualGrid_Guess = log10( abs( PB_Grid_Num(rGrid, ParaGrid, ParaSystem) ) );

%         Elevation = wrapToPi( ParaGrid.ThetaGrid + pi/2 );
        Elevation = pi/2 - ParaGrid.ThetaGrid;
        [XNoseMB, YNoseMB, ZNoseMB] = sph2cart(ParaGrid.PhiGrid, Elevation, rGuess .* ParaSystem.r0 ./ ParaSystem.Rp);
        YNose = XNoseMB;
        ZNose = YNoseMB;
        XNose = ZNoseMB;

        exNose = ParaSystem.NoseDirection;

        MTilted = ParaSystem.Tilt.RotationPhi*( ParaSystem.Tilt.RotationTheta* [0; 0; -1]);
        ezNose_temp = -(MTilted - dot(exNose, MTilted)*exNose);
        ezNose = ezNose_temp ./ norm(ezNose_temp);

        eyNose = cross(ezNose, exNose);

        P_CartNose2KSM  = [exNose, eyNose, ezNose];

        X_KSM = P_CartNose2KSM(1,1) .* XNose + P_CartNose2KSM(1,2) .* YNose + P_CartNose2KSM(1,3) .* ZNose;
        Y_KSM = P_CartNose2KSM(2,1) .* XNose + P_CartNose2KSM(2,2) .* YNose + P_CartNose2KSM(2,3) .* ZNose;
        Z_KSM = P_CartNose2KSM(3,1) .* XNose + P_CartNose2KSM(3,2) .* YNose + P_CartNose2KSM(3,3) .* ZNose;

        figure;
        hold on
        Surface = surf(X_KSM, Y_KSM, Z_KSM, X_KSM);
        Contours = plot3(X_KSM, Y_KSM, Z_KSM);

        set(Surface, 'FaceAlpha', 0.6, 'FaceColor', 'texturemap', 'Cdata', ResidualGrid_Guess);
        set(Surface, 'FaceAlpha', 0.5, 'EdgeAlpha', 0)
        set(Contours, 'LineWidth', 1, 'Color', [0 0 0 0.6])
        cb = colorbar( 'XTickLabel',{'0.1','1','10','100'}, ...
                       'XTick', [-3, -2, -1, 0] );    

        YNose = -YNose;           
        X_KSM = P_CartNose2KSM(1,1) .* XNose + P_CartNose2KSM(1,2) .* YNose + P_CartNose2KSM(1,3) .* ZNose;
        Y_KSM = P_CartNose2KSM(2,1) .* XNose + P_CartNose2KSM(2,2) .* YNose + P_CartNose2KSM(2,3) .* ZNose;
        Z_KSM = P_CartNose2KSM(3,1) .* XNose + P_CartNose2KSM(3,2) .* YNose + P_CartNose2KSM(3,3) .* ZNose;
        Surface2 = surf(X_KSM, Y_KSM, Z_KSM, X_KSM);
        Contours2 = plot3(X_KSM, Y_KSM, Z_KSM);

        set(Surface2, 'FaceAlpha', 0.6, 'FaceColor', 'texturemap', 'Cdata', ResidualGrid_Guess);
        cb = colorbar( 'XTickLabel',{'0.1','1','10','100'}, ...
                       'XTick', [-3, -2, -1, 0] );

        caxis([-3 0])


        set(Surface2, 'FaceAlpha', 0.5, 'EdgeAlpha', 0)
        set(Contours2, 'LineWidth', 1, 'Color', [0 0 0 0.6])
        axis equal
        axis([-30 30 -30 30 -30 30])

%         SurfaceTemplate.rInterp = griddedInterpolant(ParaGrid.PhiGrid, ParaGrid.ThetaGrid, rGuess);


%% TEST TEMPLATE
load('SurfaceTemplate.mat', 'SurfaceTemplate')

        ParaGrid.DeltaPhi = 5 .* DegToRad;
        ParaGrid.DeltaTheta = 5 .* DegToRad;
        ParaGrid.ThetaList = (0: ParaGrid.DeltaTheta * RadToDeg :90).' .* DegToRad;
        ParaGrid.PhiList = (-90: ParaGrid.DeltaPhi * RadToDeg : 90).' .* DegToRad;
        [ParaGrid.PhiGrid, ParaGrid.ThetaGrid] = ndgrid(ParaGrid.PhiList, ParaGrid.ThetaList);
        
rGrid = (norm(ParaSystem.Nose.KSMU_Rp) .* ParaSystem.Rp ./ ParaSystem.r0) .* SurfaceTemplate.rInterp(ParaGrid.PhiGrid, ParaGrid.ThetaGrid);
rGuess = rGrid;

% Elevation = wrapToPi( ParaGrid.ThetaGrid + pi/2 );
Elevation = ( pi/2 - ParaGrid.ThetaGrid );
[Y, Z, X] = sph2cart(ParaGrid.PhiGrid, Elevation, rGrid .* ParaSystem.r0 ./ ParaSystem.Rp);

BCAN_X = BCAN_Rotated.XInterp( X, Y, Z );
BCAN_Y = BCAN_Rotated.YInterp( X, Y, Z );
BCAN_Z = BCAN_Rotated.ZInterp( X, Y, Z );

% figure;
%     pcolor((BCAN_Z))
%     colorbar
%     
% figure; 
%     quiver3(X, Y, Z, BCAN_X, BCAN_Y, BCAN_Z)
%     axis equal

%% Computing CAN Disk X-Y-Z KSM Components
%{
mu0I = (-60.4*10^(-9) / ParaSystem.b0);
a = 8;
b = 15.5;
D = 3;
CAN_DiskParameters = [ mu0I a b D];

rList_Rp = (10 : 1 : 35) ;
PhiList = (-90:3:90) .' .* DegToRad;
ThetaList = (0:3:90) .' .* DegToRad;

[GridPhi, GridTheta, GridR] = ndgrid( PhiList, ThetaList, rList_Rp );

BCAN = B_CAN_Corrected(CAN_DiskParameters, GridR, GridTheta, GridPhi);
BCAN.XInterp = griddedInterpolant(GridPhi, GridTheta, GridR, BCAN.X);
BCAN.YInterp = griddedInterpolant(GridPhi, GridTheta, GridR, BCAN.Y);
BCAN.ZInterp = griddedInterpolant(GridPhi, GridTheta, GridR, BCAN.Z);

ParaSystem.BCAN = BCAN;
%}


%{
XGrid = GridR .* cos(GridTheta);
YGrid = GridR .* sin(GridTheta) .* cos(GridPhi);
ZGrid = GridR .* sin(GridTheta) .* sin(GridPhi);

figure;
hold on
    PlotQuiver = quiver3(XGrid, YGrid, ZGrid, BCAN_XInterp(GridPhi, GridTheta, GridR), BCAN_YInterp(GridPhi, GridTheta, GridR), BCAN_ZInterp(GridPhi, GridTheta, GridR));
    axis equal
    set(PlotQuiver, 'AutoScaleFactor', 4, 'LineWidth', 1, 'MaxHeadSize', 1);
    
figure;
hold on
BCAN_r = sqrt(BCAN_XInterp(GridPhi, GridTheta, GridR).^2 + BCAN_YInterp(GridPhi, GridTheta, GridR).^2 + BCAN_ZInterp(GridPhi, GridTheta, GridR).^2);
    scatter3( XGrid(:), ZGrid(:), BCAN_r(:)  )
%}


    %% -------------------------------------------------------------------------------
    % Optimising Towards Pressure Balance: Iterative Correction

        MaxIter = 3;
        MaxEval = 5e5;
        ParaGrid.epsilon = 1e-10;
        [PBGrid_Func] = @(r) PB_Grid(r, ParaGrid, ParaSystem);
        Algorithms = ["levenberg-marquardt"; "Trust-Region-Dogleg"; "Trust-Region" ];
        opts = optimoptions('fsolve', 'Display','iter', 'Maxiter', MaxIter, 'MaxFunctionEvaluations', MaxEval,...
                            'SpecifyObjectiveGradient', false, 'CheckGradients', false, 'FiniteDifferenceType', 'central', 'Algorithm', Algorithms(1),...
                            'ScaleProblem', 'none');

        rGuess = rGrid;
        rSol_Vec = fsolve(PBGrid_Func, rGuess, opts);
        rSol = reshape(rSol_Vec, [length(ParaGrid.PhiList), length(ParaGrid.ThetaList)]);
        rGrid = rSol;

    %     [Grid, test] = PBGrid_Func(rGuess)


    rInterp = griddedInterpolant(ParaGrid.PhiGrid, ParaGrid.ThetaGrid, rSol);
    rGuess = rInterp(ParaGrid.PhiGrid, ParaGrid.ThetaGrid);
    
    
    
%     rSol(end-1, :) = (1/2) .* ( rSol(end, :) +  rSol(end-2, :) );
    
    
    %{
%% Introduction BShielding

% Bx Interp
    BShielding = ShieldingField(rGrid, ParaGrid, ParaSystem);
    rList = rGrid(:);
        
    BXGrid = BShielding.XInterp(ParaGrid.PhiGrid, ParaGrid.ThetaGrid);
    BXList = BXGrid(:);
    
    [~,rListsort]=sort(rList);
    rListSorted  = rList(rListsort) ;
    BXListSorted  = BXList(rListsort) ;

    repeats = diff(rListSorted) == 0; 
    rListSorted(repeats) = [];
    BXListSorted(repeats) = [];

    BShielding.XInterp_r = griddedInterpolant(rListSorted, BXListSorted);
    
% By Interp
    BYGrid = BShielding.YInterp(ParaGrid.PhiGrid, ParaGrid.ThetaGrid);
    BYList = BYGrid(:);
    BYListSorted  = BYList(rListsort) ;
    BYListSorted(repeats) = [];
    
    BShielding.YInterp_r = griddedInterpolant(rListSorted, BYListSorted);
        
    
% Bz Interp
    BZGrid = BShielding.ZInterp(ParaGrid.PhiGrid, ParaGrid.ThetaGrid);
    BZList = BZGrid(:);
    BZListSorted  = BZList(rListsort) ;
    BZListSorted(repeats) = [];
    
    BShielding.ZInterp_r = griddedInterpolant(rListSorted, BZListSorted);

ParaSystem.BShielding = BShielding;

%}






%     griddedInterpolant(rGrid(:), BShielding.X(:))
% 
% DeltaR =  1 .* ParaSystem.Rp ./ ParaSystem.r0;
% rList = [min(rSol(:)) : DeltaR : max(rSol(:))] .';
% [GridPhi, GridTheta, GridR] = ndgrid(ParaGrid.PhiList, ParaGrid.ThetaList, rList);
% Grid_Bx = 0.*GridR;
% Grid_By = 0.*GridR;
% Grid_Bz = 0.*GridR;
% 
% 
% for kr = 1:length(rList)
%     rGrid = GridR(:, :, kr);
%     BShielding = ShieldingField(rGrid, ParaGrid, ParaSystem);
%     Grid_Bx(:, :, kr) = BShielding.XInterp(GridPhi(:, :, kr), GridTheta(:, :, kr));
%     Grid_By(:, :, kr) = BShielding.YInterp(GridPhi(:, :, kr), GridTheta(:, :, kr));
%     Grid_Bz(:, :, kr) = BShielding.ZInterp(GridPhi(:, :, kr), GridTheta(:, :, kr));
% end
% 
% BShielding.XInterp_3 = griddedInterpolant(GridPhi, GridTheta, GridR, Grid_Bx);
% BShielding.YInterp_3 = griddedInterpolant(GridPhi, GridTheta, GridR, Grid_By);
% BShielding.ZInterp_3 = griddedInterpolant(GridPhi, GridTheta, GridR, Grid_Bz);
% 
% ParaSystem.BShielding = BShielding;
% 
% 
% DeltaR =  1 .* ParaSystem.Rp ./ ParaSystem.r0;
% rList = [min(rSol(:)) : DeltaR : max(rSol(:))] .';
% 
%      BShielding = ShieldingField(rList, ParaGrid, ParaSystem);


%% TEST
% figure
%     rSol = rGrid;
%     method = 'gaussian';
%     window1 = 2;
% 	window2 = 2;
% 
%         rSol = smoothdata(rSol, 2, method, window1);
%         rSol = smoothdata(rSol, 1, method, window2);
% 

           
    %% TEST: Comparing Jacobian
   
    ResidualGrid_Sol = log10( abs( PB_Grid(rSol, ParaGrid, ParaSystem) ) );
%     figure;
%     pcolor(ResidualGrid_Sol)
%     cb = colorbar( 'XTickLabel',{'0.1','1','10','100'}, ...
%                    'XTick', [-3, -2, -1, 0] );
%     caxis([-3 0])
% %     

%% -------------------------------------------------------------------------------
% Optimising Towards Pressure Balance: Conversion of Surface into KSM Frame

%         Elevation = wrapToPi( ParaGrid.ThetaGrid + pi/2 );
        Elevation = pi/2 - ParaGrid.ThetaGrid;
        [XNoseMB, YNoseMB, ZNoseMB] = sph2cart(ParaGrid.PhiGrid, Elevation, rSol .* ParaSystem.r0 ./ ParaSystem.Rp);
%     [XNoseMB, YNoseMB, ZNoseMB] = sph2cart(ParaGrid.PhiGrid, pi/2-ParaGrid.ThetaGrid, rSol .* ParaSystem.r0 ./ ParaSystem.Rp);
    YNose = XNoseMB;
    ZNose = YNoseMB;
    XNose = ZNoseMB;

    exNose = ParaSystem.NoseDirection;

    MTilted = ParaSystem.Tilt.RotationPhi*( ParaSystem.Tilt.RotationTheta* [0; 0; -1]);
    ezNose_temp = -(MTilted - dot(exNose, MTilted)*exNose);
    ezNose = ezNose_temp ./ norm(ezNose_temp);

    eyNose = cross(ezNose, exNose);

    P_CartNose2KSM  = [exNose, eyNose, ezNose];

    X_KSM = P_CartNose2KSM(1,1) .* XNose + P_CartNose2KSM(1,2) .* YNose + P_CartNose2KSM(1,3) .* ZNose;
    Y_KSM = P_CartNose2KSM(2,1) .* XNose + P_CartNose2KSM(2,2) .* YNose + P_CartNose2KSM(2,3) .* ZNose;
    Z_KSM = P_CartNose2KSM(3,1) .* XNose + P_CartNose2KSM(3,2) .* YNose + P_CartNose2KSM(3,3) .* ZNose;
    
    figure;
    hold on
    Surface = surf(X_KSM, Y_KSM, Z_KSM, X_KSM);
    XValue = 16;
    Contours = plot3(X_KSM(XValue-0.5<XNose<XValue+0.5 & abs(ZNose) > 1), Y_KSM(XValue-0.5<XNose<XValue+0.5 & abs(ZNose) > 1), Z_KSM(XValue-0.5<XNose<XValue+0.5 & abs(ZNose)>1));

    set(Surface, 'FaceAlpha', 0.6, 'FaceColor', 'texturemap', 'Cdata', ResidualGrid_Sol);
    cb = colorbar( 'XTickLabel',{'0.1','1','10','100'}, ...
                   'XTick', [-3, -2, -1, 0] );
    caxis([-3 0])


    set(Surface, 'FaceAlpha', 0.5, 'EdgeAlpha', 0)
    set(Contours, 'LineWidth', 1, 'Color', [0 0 0 0.6])
    axis equal
    
    
    YNose = -YNose;           
    X_KSM = P_CartNose2KSM(1,1) .* XNose + P_CartNose2KSM(1,2) .* YNose + P_CartNose2KSM(1,3) .* ZNose;
    Y_KSM = P_CartNose2KSM(2,1) .* XNose + P_CartNose2KSM(2,2) .* YNose + P_CartNose2KSM(2,3) .* ZNose;
    Z_KSM = P_CartNose2KSM(3,1) .* XNose + P_CartNose2KSM(3,2) .* YNose + P_CartNose2KSM(3,3) .* ZNose;
    Surface2 = surf(X_KSM, Y_KSM, Z_KSM, X_KSM);
    Contours2 = plot3(X_KSM, Y_KSM, Z_KSM);
    
    set(Surface2, 'FaceAlpha', 0.6, 'FaceColor', 'texturemap', 'Cdata', ResidualGrid_Sol);
    cb = colorbar( 'XTickLabel',{'0.1','1','10','100'}, ...
                   'XTick', [-3, -2, -1, 0] );
               
    caxis([-3 0])


    set(Surface2, 'FaceAlpha', 0.5, 'EdgeAlpha', 0)
    set(Contours2, 'LineWidth', 1, 'Color', [0 0 0 0.6])
    axis equal
    axis([-30 30 -30 30 -30 30])
    
    
    
%  [X, Y, Z] =    surfnorm(X_KSM, Y_KSM, Z_KSM)
    
%     
%     figure;
%     hold on
%     XNose_Interp = griddedInterpolant(ParaGrid.PhiGrid, ParaGrid.ThetaGrid, XNose);
%     YNose_Interp = griddedInterpolant(ParaGrid.PhiGrid, ParaGrid.ThetaGrid, YNose);
%     ZNose_Interp = griddedInterpolant(ParaGrid.PhiGrid, ParaGrid.ThetaGrid, ZNose);
%     
%     Phi = (-90:3:90) .* pi/180;
%     Theta = (0:1:100) .* pi/180;
%     [PhiGrid, ThetaGrid] = ndgrid(Phi, Theta);
%     
%     XNose  = XNose_Interp(PhiGrid, ThetaGrid);
%     YNose  = YNose_Interp(PhiGrid, ThetaGrid);
%     ZNose  = ZNose_Interp(PhiGrid, ThetaGrid);
% 
%     X_KSM = P_CartNose2KSM(1,1) .* XNose + P_CartNose2KSM(1,2) .* YNose + P_CartNose2KSM(1,3) .* ZNose;
%     Y_KSM = P_CartNose2KSM(2,1) .* XNose + P_CartNose2KSM(2,2) .* YNose + P_CartNose2KSM(2,3) .* ZNose;
%     Z_KSM = P_CartNose2KSM(3,1) .* XNose + P_CartNose2KSM(3,2) .* YNose + P_CartNose2KSM(3,3) .* ZNose;
%             
%     surf( X_KSM, Y_KSM, Z_KSM)
%     
%     XNose  = XNose_Interp(PhiGrid, ThetaGrid);
%     YNose  = -YNose_Interp(PhiGrid, ThetaGrid);
%     ZNose  = ZNose_Interp(PhiGrid, ThetaGrid);
% 
%     X_KSM = P_CartNose2KSM(1,1) .* XNose + P_CartNose2KSM(1,2) .* YNose + P_CartNose2KSM(1,3) .* ZNose;
%     Y_KSM = P_CartNose2KSM(2,1) .* XNose + P_CartNose2KSM(2,2) .* YNose + P_CartNose2KSM(2,3) .* ZNose;
%     Z_KSM = P_CartNose2KSM(3,1) .* XNose + P_CartNose2KSM(3,2) .* YNose + P_CartNose2KSM(3,3) .* ZNose;
%             
%     surf( X_KSM, Y_KSM, Z_KSM)
%     
%     axis equal
%     
    
  %% TEST: Nose Determination
  
ThetaTilt = 27.* DegToRad;
PhiList = (0: 1 : 0 + 360).' .* DegToRad;
% PhiList = (-180: 1 : 0 + 180).' .* DegToRad;
ThetaList = 0.*PhiList + ThetaTilt;
ThetaSolList = 0.*PhiList;

rList = 0 .* PhiList + 1;


for kPhi = 1:length(PhiList)
    Phi = PhiList(kPhi);
    Theta = ThetaList(1);
        Func = @(theta) sec(theta).^3+(-1).*(sin(Phi).^2.*sin(Theta).^2+(cos(Theta).*((-1) ...
  +3.*sin(theta).^2)+3.*cos(Phi).*cos(theta).*sin(theta).*sin(Theta) ...
  ).^2).^(1/2);

    if kPhi == 1
        ThetaGuess = ThetaTilt;
    else
        ThetaGuess = ThetaSolList(kPhi-1);
    end
    
    thetaSol = fsolve(Func, ThetaGuess);
    ThetaSolList(kPhi) = thetaSol;
    
        ParaSystem.Tilt.Phi = Phi;
        ParaSystem.Tilt.RotationPhi = [ cos(ParaSystem.Tilt.Phi) -sin(ParaSystem.Tilt.Phi) 0 ;...
                                    sin(ParaSystem.Tilt.Phi)  cos(ParaSystem.Tilt.Phi) 0 ;...
                                    0              0             1];

    PBNose = @(r, NoseTheta) PB_Nose(r, NoseTheta, ParaSystem);
    rSol = fsolve( @(r) PBNose(r, thetaSol), 1);
    rList(kPhi) = rSol;

end

PhiNMM = acos(sin(PhiList) .* sin(ThetaList));

rList_Rp = rList .* ParaSystem.r0 ./ ParaSystem.Rp;
XList = rList_Rp .* cos(ThetaSolList);
YList = rList_Rp .* sin(ThetaSolList) .* cos(PhiNMM);
ZList = rList_Rp .* sin(ThetaSolList) .* sin(PhiNMM);

%     XList = P_CartNose2KSM(1,1) .* XList + P_CartNose2KSM(1,2) .* YList + P_CartNose2KSM(1,3) .* ZList;
%     YList = P_CartNose2KSM(2,1) .* XList + P_CartNose2KSM(2,2) .* YList + P_CartNose2KSM(2,3) .* ZList;
%     ZList = P_CartNose2KSM(3,1) .* XList + P_CartNose2KSM(3,2) .* YList + P_CartNose2KSM(3,3) .* ZList;

figure;
    hold on
    % Nose locus
        NoseLocus = scatter3(XList(2:end), YList(2:end), ZList(2:end), 40, PhiList(2:end) * RadToDeg, 'filled');
        Color_NMMPlane = cbrewer2('Greens', 40);
        Opacity = 0.5;
%         set(NoseLocus, 'MarkerFaceColor', Color_NMMPlane(10,:), 'MarkerFaceAlpha', Opacity)
        colorbar
        
    C = 5;
    axis([15, 20, -C, C, -C, C])
    axis equal
    ColorMap = vertcat( cbrewer2( 'BuPu', 100 ), flipud(cbrewer2( 'BuPu', 100 )) );
    colormap(ColorMap)    
    
        % ---- Plot options
        axes = gca;
        set(axes, 'FontName'   , 'Helvetica',...
            'XAxisLocation', 'top',...
            'YAxisLocation', 'left',...
            'FontSize', 16,...
            'LineWidth'   , 1 ...
            );
        Ylabel = ylabel('X (R_p)');
        Xlabel = xlabel('Y (R_p)');
        
        Zlabel = zlabel('Z (R_p)');
        
        
    
    %{
    % Surfaces
    X_KSM = P_CartNose2KSM(1,1) .* XNose + P_CartNose2KSM(1,2) .* YNose + P_CartNose2KSM(1,3) .* ZNose;
    Y_KSM = P_CartNose2KSM(2,1) .* XNose + P_CartNose2KSM(2,2) .* YNose + P_CartNose2KSM(2,3) .* ZNose;
    Z_KSM = P_CartNose2KSM(3,1) .* XNose + P_CartNose2KSM(3,2) .* YNose + P_CartNose2KSM(3,3) .* ZNose;
    Surface = surf(X_KSM, Y_KSM, Z_KSM, X_KSM);
    Contours = plot3(X_KSM, Y_KSM, Z_KSM);
    
    YNose = -YNose;           
    X_KSM = P_CartNose2KSM(1,1) .* XNose + P_CartNose2KSM(1,2) .* YNose + P_CartNose2KSM(1,3) .* ZNose;
    Y_KSM = P_CartNose2KSM(2,1) .* XNose + P_CartNose2KSM(2,2) .* YNose + P_CartNose2KSM(2,3) .* ZNose;
    Z_KSM = P_CartNose2KSM(3,1) .* XNose + P_CartNose2KSM(3,2) .* YNose + P_CartNose2KSM(3,3) .* ZNose;
    Surface2 = surf(X_KSM, Y_KSM, Z_KSM, X_KSM);
    Contours2 = plot3(X_KSM, Y_KSM, Z_KSM);

    set([ Contours, Contours2], 'LineWidth', 1, 'Color', [0 0 0 0.4])
    set([Surface, Surface2], 'FaceColor', 'black', 'FaceAlpha', 0.3, 'EdgeAlpha', 0)
    ColorMap = vertcat( cbrewer2( 'BuPu', 100 ), flipud(cbrewer2( 'BuPu', 100 )) );
    colormap(ColorMap)    
    
    C = 20;
    axis([-C, C, -C, C, -C, C])
    %}
    
% figure;
% scatter( PhiList.* RadToDeg, ThetaSolList.* RadToDeg)
% axis equal



%% TEST: Interpolate MP Surfaces

PhiTiltList = [0:10:180];
[GridPhi, GridTheta, GridPhiTilt] = ndgrid(ParaGrid.PhiList, ParaGrid.ThetaList, PhiTiltList);

GridR = 0*GridPhi;

for kPhi = 0:1:18
    Phi = 10.*kPhi;
    FileName = horzcat( 'ResultSurfaces/Phi', num2str(Phi), '.mat');
    load('ResultSurfaces/Phi30.mat', 'rInterp', 'rSol');
    GridR(:, :, kPhi+1) = rSol;
    
end
rInterp = griddedInterpolant(GridPhi, GridTheta, GridPhiTilt, GridR, 'spline');


figure;
hold on
for PhiTilt = 0:2:180
% for PhiTilt = 0

    cla

PhiTilt
[GridPhi, GridTheta, GridPhiTilt] = ndgrid(ParaGrid.PhiList, ParaGrid.ThetaList, PhiTilt);
rSurface = rInterp(GridPhi, GridTheta, GridPhiTilt);


% Building Surface


    ParaSystem.Tilt.Phi = PhiTilt * DegToRad;
    ParaSystem.Tilt.RotationPhi = [ cos(ParaSystem.Tilt.Phi) -sin(ParaSystem.Tilt.Phi) 0 ;...
                                sin(ParaSystem.Tilt.Phi)  cos(ParaSystem.Tilt.Phi) 0 ;...
                                0              0             1];
                        
    AlphaDipole = pi/2 - acos(sin(ParaSystem.Tilt.Phi)*sin(ParaSystem.Tilt.Theta));
    ParaSystem.AlphaDipole = AlphaDipole;
    ParaSystem.Tilt.RotationAlphaDipole = [ 1 0                 0                 ;...
                                            0 cos(AlphaDipole) -sin(AlphaDipole)  ;...
                                            0 sin(AlphaDipole)  cos(AlphaDipole)  ];

    Plot3D_ini(Rotation, ParaSystem);
% 	view(110, -10)
	view(90, 0)

                            
    [XNoseMB, YNoseMB, ZNoseMB] = sph2cart(ParaGrid.PhiGrid, pi/2-ParaGrid.ThetaGrid, rSurface .* ParaSystem.r0 ./ ParaSystem.Rp);
    YNose = XNoseMB;
    ZNose = YNoseMB;
    XNose = ZNoseMB;

    exNose = ParaSystem.NoseDirection;

    MTilted = ParaSystem.Tilt.RotationPhi*( ParaSystem.Tilt.RotationTheta* [0; 0; -1]);
    ezNose_temp = -(MTilted - dot(exNose, MTilted)*exNose);
    ezNose = ezNose_temp ./ norm(ezNose_temp);

    eyNose = cross(ezNose, exNose);

    P_CartNose2KSM  = [exNose, eyNose, ezNose];

    X_KSM = P_CartNose2KSM(1,1) .* XNose + P_CartNose2KSM(1,2) .* YNose + P_CartNose2KSM(1,3) .* ZNose;
    Y_KSM = P_CartNose2KSM(2,1) .* XNose + P_CartNose2KSM(2,2) .* YNose + P_CartNose2KSM(2,3) .* ZNose;
    Z_KSM = P_CartNose2KSM(3,1) .* XNose + P_CartNose2KSM(3,2) .* YNose + P_CartNose2KSM(3,3) .* ZNose;
    
%     figure;
    hold on
    Surface = surf(X_KSM, Y_KSM, Z_KSM, X_KSM);
    XValue = 16;
    Contours = plot3(X_KSM(XValue-0.5<XNose<XValue+0.5 & abs(ZNose) > 1), Y_KSM(XValue-0.5<XNose<XValue+0.5 & abs(ZNose) > 1), Z_KSM(XValue-0.5<XNose<XValue+0.5 & abs(ZNose)>1));

    set(Surface, 'FaceAlpha', 0.6, 'FaceColor', 'texturemap', 'Cdata', ResidualGrid_Sol);
    cb = colorbar( 'XTickLabel',{'0.1','1','10','100'}, ...
                   'XTick', [-3, -2, -1, 0] );
    caxis([-3 0])


    set(Surface, 'FaceAlpha', 0.3, 'EdgeAlpha', 0)
    set(Contours, 'LineWidth', 1, 'Color', [0 0 0 0.6])
    axis equal
    
    
    YNose = -YNose;           
    X_KSM = P_CartNose2KSM(1,1) .* XNose + P_CartNose2KSM(1,2) .* YNose + P_CartNose2KSM(1,3) .* ZNose;
    Y_KSM = P_CartNose2KSM(2,1) .* XNose + P_CartNose2KSM(2,2) .* YNose + P_CartNose2KSM(2,3) .* ZNose;
    Z_KSM = P_CartNose2KSM(3,1) .* XNose + P_CartNose2KSM(3,2) .* YNose + P_CartNose2KSM(3,3) .* ZNose;
    Surface2 = surf(X_KSM, Y_KSM, Z_KSM, X_KSM);
    Contours2 = plot3(X_KSM, Y_KSM, Z_KSM);
    
    set(Surface2, 'FaceAlpha', 0.6, 'FaceColor', 'texturemap', 'Cdata', ResidualGrid_Sol);
    cb = colorbar( 'XTickLabel',{'0.1','1','10','100'}, ...
                   'XTick', [-3, -2, -1, 0] );
               
    caxis([-3 0])


    set(Surface2, 'FaceAlpha', 0.3, 'EdgeAlpha', 0)
    set(Contours2, 'LineWidth', 1, 'Color', [0 0 0 0.6])
    axis equal
    axis([-30 30 -30 30 -30 30])
    
    View =  [ [ -0.4684    0.8835         0   -0.2075 ]; ...
              [ -0.1052   -0.0558    0.9929   -0.4159 ]; ...
              [ -0.8772   -0.4651   -0.1191    9.3910  ]; ...
             [      0         0         0    1.0000  ] ];
         
         view( 110, 10)
%     camlight('left')
     pause(0.001);


end














load('ResultSurfaces/Phi0.mat')
rSol_0 = rSol;
GridTheta = ParaGrid.ThetaGrid;
GridPhi = ParaGrid.PhiGrid;

load('ResultSurfaces/Phi30.mat')
rSol_30 = rSol;

load('ResultSurfaces/Phi60.mat')
rSol_60 = rSol;

load('ResultSurfaces/Phi90.mat')
rSol_90 = rSol;

load('ResultSurfaces/Phi120.mat')
rSol_120 = rSol;



[GridPhi, GridTheta, GridPhiTilt] = ndgrid(ParaGrid.PhiList, ParaGrid.ThetaList, PhiTiltList);

GridR = 0*GridPhi;
GridR(:, :, 1) = rSol_0;
GridR(:, :, 2) = rSol_30;
GridR(:, :, 3) = rSol_60;
GridR(:, :, 4) = rSol_90;
GridR(:, :, 5) = rSol_120;


rInterp = griddedInterpolant(GridPhi, GridTheta, GridPhiTilt, GridR, 'spline');


figure;
for PhiTilt = 0:2:120

    cla

PhiTilt
[GridPhi, GridTheta, GridPhiTilt] = ndgrid(ParaGrid.PhiList, ParaGrid.ThetaList, PhiTilt);
rSurface = rInterp(GridPhi, GridTheta, GridPhiTilt);


% Building Surface


    ParaSystem.Tilt.Phi = PhiTilt * DegToRad;
    ParaSystem.Tilt.RotationPhi = [ cos(ParaSystem.Tilt.Phi) -sin(ParaSystem.Tilt.Phi) 0 ;...
                                sin(ParaSystem.Tilt.Phi)  cos(ParaSystem.Tilt.Phi) 0 ;...
                                0              0             1];
                        
    AlphaDipole = pi/2 - acos(sin(ParaSystem.Tilt.Phi)*sin(ParaSystem.Tilt.Theta));
    ParaSystem.AlphaDipole = AlphaDipole;
    ParaSystem.Tilt.RotationAlphaDipole = [ 1 0                 0                 ;...
                                            0 cos(AlphaDipole) -sin(AlphaDipole)  ;...
                                            0 sin(AlphaDipole)  cos(AlphaDipole)  ];

    Plot3D_ini(Rotation, ParaSystem);
% 	view(110, -10)
	view(90, 0)

                            
    [XNoseMB, YNoseMB, ZNoseMB] = sph2cart(ParaGrid.PhiGrid, pi/2-ParaGrid.ThetaGrid, rSurface .* ParaSystem.r0 ./ ParaSystem.Rp);
    YNose = XNoseMB;
    ZNose = YNoseMB;
    XNose = ZNoseMB;

    exNose = ParaSystem.NoseDirection;

    MTilted = ParaSystem.Tilt.RotationPhi*( ParaSystem.Tilt.RotationTheta* [0; 0; -1]);
    ezNose_temp = -(MTilted - dot(exNose, MTilted)*exNose);
    ezNose = ezNose_temp ./ norm(ezNose_temp);

    eyNose = cross(ezNose, exNose);

    P_CartNose2KSM  = [exNose, eyNose, ezNose];

    X_KSM = P_CartNose2KSM(1,1) .* XNose + P_CartNose2KSM(1,2) .* YNose + P_CartNose2KSM(1,3) .* ZNose;
    Y_KSM = P_CartNose2KSM(2,1) .* XNose + P_CartNose2KSM(2,2) .* YNose + P_CartNose2KSM(2,3) .* ZNose;
    Z_KSM = P_CartNose2KSM(3,1) .* XNose + P_CartNose2KSM(3,2) .* YNose + P_CartNose2KSM(3,3) .* ZNose;
    
%     figure;
    hold on
    Surface = surf(X_KSM, Y_KSM, Z_KSM, X_KSM);
    XValue = 16;
    Contours = plot3(X_KSM(XValue-0.5<XNose<XValue+0.5 & abs(ZNose) > 1), Y_KSM(XValue-0.5<XNose<XValue+0.5 & abs(ZNose) > 1), Z_KSM(XValue-0.5<XNose<XValue+0.5 & abs(ZNose)>1));

    set(Surface, 'FaceAlpha', 0.6, 'FaceColor', 'texturemap', 'Cdata', ResidualGrid_Sol);
    cb = colorbar( 'XTickLabel',{'0.1','1','10','100'}, ...
                   'XTick', [-3, -2, -1, 0] );
    caxis([-3 0])


    set(Surface, 'FaceAlpha', 0.3, 'EdgeAlpha', 0)
    set(Contours, 'LineWidth', 1, 'Color', [0 0 0 0.6])
    axis equal
    
    
    YNose = -YNose;           
    X_KSM = P_CartNose2KSM(1,1) .* XNose + P_CartNose2KSM(1,2) .* YNose + P_CartNose2KSM(1,3) .* ZNose;
    Y_KSM = P_CartNose2KSM(2,1) .* XNose + P_CartNose2KSM(2,2) .* YNose + P_CartNose2KSM(2,3) .* ZNose;
    Z_KSM = P_CartNose2KSM(3,1) .* XNose + P_CartNose2KSM(3,2) .* YNose + P_CartNose2KSM(3,3) .* ZNose;
    Surface2 = surf(X_KSM, Y_KSM, Z_KSM, X_KSM);
    Contours2 = plot3(X_KSM, Y_KSM, Z_KSM);
    
    set(Surface2, 'FaceAlpha', 0.6, 'FaceColor', 'texturemap', 'Cdata', ResidualGrid_Sol);
    cb = colorbar( 'XTickLabel',{'0.1','1','10','100'}, ...
                   'XTick', [-3, -2, -1, 0] );
               
    caxis([-3 0])


    set(Surface2, 'FaceAlpha', 0.3, 'EdgeAlpha', 0)
    set(Contours2, 'LineWidth', 1, 'Color', [0 0 0 0.6])
    axis equal
    axis([-30 30 -30 30 -30 30])
    
     pause(0.005);


end
