
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
        Conditions_PhiTilt_rMP_Psw = [ 0, 25, NaN ];
        ParaSystem = System_ini(Planet, Inclination, Conditions_PhiTilt_rMP_Psw);
        ParaSystem.BShielding.XVal = 0;
        ParaSystem.BShielding.YVal = 0;
        ParaSystem.BShielding.ZVal = 0;


    %% ------------------------------------------------------------------------
    % Computing ROTATED CAN Disk X-Y-Z KSM Components FROM KSM coordinates

    % ----  Setting up Cartesian grid for CAN Disk Field Interpolation
        ParaSystem.OnOff_CAN = 1;
        CAN_DiskParameters = ParaSystem.CAN_DiskParameters;

        XList_KSM = (-36 : 2 : 36) .';   
        YList_KSM = (-36 : 2 : 36) .';   
        ZList_KSM = (-36 : 2 : 36) .';   
        
        XList_KSM( abs(XList_KSM) <= 5 ) = [];        
        YList_KSM( abs(YList_KSM) <= 5 ) = [];        
        ZList_KSM( abs(ZList_KSM) <= 5 ) = [];        

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

    % ----  Rotating to initial tilted-dipole case: Position_Tilted = [XKSM; YKSM; ZKSM]
        BCAN_Rotated.X = P_AlignedToRot(1,1) .* BCAN_Aligned.X + P_AlignedToRot(1,2) .* BCAN_Aligned.Y  + P_AlignedToRot(1,3) .* BCAN_Aligned.Z;
        BCAN_Rotated.Y = P_AlignedToRot(2,1) .* BCAN_Aligned.X + P_AlignedToRot(2,2) .* BCAN_Aligned.Y  + P_AlignedToRot(2,3) .* BCAN_Aligned.Z;
        BCAN_Rotated.Z = P_AlignedToRot(3,1) .* BCAN_Aligned.X + P_AlignedToRot(3,2) .* BCAN_Aligned.Y  + P_AlignedToRot(3,3) .* BCAN_Aligned.Z;

        BCAN_Rotated.XInterp = griddedInterpolant(XGrid_KSM, YGrid_KSM, ZGrid_KSM, BCAN_Rotated.X, 'linear');
        BCAN_Rotated.YInterp = griddedInterpolant(XGrid_KSM, YGrid_KSM, ZGrid_KSM, BCAN_Rotated.Y, 'linear');
        BCAN_Rotated.ZInterp = griddedInterpolant(XGrid_KSM, YGrid_KSM, ZGrid_KSM, BCAN_Rotated.Z, 'linear');
        
    % ----  Dipole Field
    
        Grid.X = XGrid_KSM;
        Grid.Y = YGrid_KSM;
        Grid.Z = ZGrid_KSM;
        BDipole = DipoleField(Grid, ParaSystem);
    
    %{
        MTilted = ParaSystem.MTilted;
        MGridX = 0*XGrid_KSM(:) + MTilted(1);
        MGridY = 0*YGrid_KSM(:) + MTilted(2);
        MGridZ = 0*ZGrid_KSM(:) + MTilted(3);
        
        Position = [ XGrid_KSM(:), YGrid_KSM(:), ZGrid_KSM(:) ] .* ParaSystem.Rp ./ ParaSystem.r0;
        Norm = sqrt (Position(:,1).^2 + Position(:,2).^2 + Position(:,3).^2 );
        er = Position ./ Norm;
        XNorm = er(:,1);
        YNorm = er(:,2);
        ZNorm = er(:,3);
        
        MdotEr = MTilted(1) .* XNorm + MTilted(2) .* YNorm + MTilted(3) .* ZNorm  ;
    
        BDipole.X = (1./Norm.^3) .* ( 3*MdotEr .* XNorm - MGridX );
        BDipole.Y = (1./Norm.^3) .* ( 3*MdotEr .* YNorm - MGridY );
        BDipole.Z = (1./Norm.^3) .* ( 3*MdotEr .* ZNorm - MGridZ );
        
        BDipole.XInterp = griddedInterpolant(XGrid_KSM, YGrid_KSM, ZGrid_KSM, reshape(BDipole.X, size(XGrid_KSM)),  'linear');
        BDipole.YInterp = griddedInterpolant(XGrid_KSM, YGrid_KSM, ZGrid_KSM, reshape(BDipole.Y, size(XGrid_KSM)),  'linear');
        BDipole.ZInterp = griddedInterpolant(XGrid_KSM, YGrid_KSM, ZGrid_KSM, reshape(BDipole.Z, size(XGrid_KSM)),  'linear');
     %}

    % ----  Optional: Visualisation of CAN-disk field 
    %{
        XList = (-20:1:20);
        YList = (-20:1:20);
        ZList = (-20:1:20);

        [X, Y, Z] = ndgrid( XList, YList, ZList );
        Bx = BDipole.XInterp(X, Y, Z);
        By = BDipole.YInterp(X, Y, Z)
        Bz = BDipole.ZInterp(X, Y, Z)
        BNorm = sqrt(Bx.^2 + By.^2 + Bz.^2)
        figure;
            PlotQuiver = quiver3(X, Y, Z, Bx./BNorm, By./BNorm, Bz./BNorm);
            set(PlotQuiver, 'AutoScaleFactor', 1, 'LineWidth', 1);
            axis equal
            axis([-20 20 0 1/2 -20 20])
    %}

    % ----  Packaging Output
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

        norm(ParaSystem.Nose.KSMU_Rp)
        scatter3(ParaSystem.Nose.KSMU_Rp(1), ParaSystem.Nose.KSMU_Rp(2), ParaSystem.Nose.KSMU_Rp(3), 100, 'filled')
        [ParaSystem.Nose.KSMU_Rp(1), ParaSystem.Nose.KSMU_Rp(2), ParaSystem.Nose.KSMU_Rp(3)]
%     close all
     
%% -------------------------------------------------------------------------------
% Building Guess-Surface

% ---- Loading previously-computed Guess Surface: Aligned Dipole, Only Dipole Field
    load('SurfaceTemplate.mat', 'SurfaceTemplate')
%     load('Phi90.mat', 'rInterp', 'rGrid')

% ----  Setting up grid for Guess Surface to optimise
    ParaGrid.DeltaPhi = 5 .* DegToRad;
    ParaGrid.DeltaTheta = 5 .* DegToRad;
    ParaGrid.ThetaList = (0: ParaGrid.DeltaTheta * RadToDeg :90).' .* DegToRad;
    ParaGrid.PhiList = (-90: ParaGrid.DeltaPhi * RadToDeg : 90).' .* DegToRad;
    [ParaGrid.PhiGrid, ParaGrid.ThetaGrid] = ndgrid(ParaGrid.PhiList, ParaGrid.ThetaList);
        
% ----  Normalising to current Nose dimensions
    rGrid = (norm(ParaSystem.Nose.KSMU_Rp) .* ParaSystem.Rp ./ ParaSystem.r0) .* SurfaceTemplate.rInterp(ParaGrid.PhiGrid, ParaGrid.ThetaGrid);
%     rGrid = (norm(ParaSystem.Nose.KSMU_Rp) .* ParaSystem.Rp ./ ParaSystem.r0) .* rInterp(ParaGrid.PhiGrid, ParaGrid.ThetaGrid);
    rGuess = rGrid;

    
    
%% -------------------------------------------------------------------------------
% Optimising Guess-Surface Towards Pressure Balance
    
% ----  Turning on/off Shielding Field
    ParaSystem.OnOff_ShieldingField = 0;

% ----  Setting up optimisation method
    MaxIter = 5;
    MaxEval = 5e5;
    [PBGrid_Func] = @(r) PB_Grid(r, ParaGrid, ParaSystem);
    Algorithms = ["levenberg-marquardt"; "Trust-Region-Dogleg"; "Trust-Region" ];
    opts = optimoptions('fsolve', 'Display','iter', 'Maxiter', MaxIter, 'MaxFunctionEvaluations', MaxEval,...
        'SpecifyObjectiveGradient', false, 'CheckGradients', false, 'FiniteDifferenceType', 'central', 'Algorithm', Algorithms(1),...
        'ScaleProblem', 'none');
    
% ----  Solving and reshaping solution
    rSol_Vec = fsolve(PBGrid_Func, rGuess, opts);
    rSol = reshape(rSol_Vec, [length(ParaGrid.PhiList), length(ParaGrid.ThetaList)]);
    rGrid = rSol;
    
% ----  Define interpolant with solution
    rInterp = griddedInterpolant(ParaGrid.PhiGrid, ParaGrid.ThetaGrid, rSol, 'spline');
    rGuess = rInterp(ParaGrid.PhiGrid, ParaGrid.ThetaGrid);
    
%     if  ParaSystem.OnOff_ShieldingField == 1
%         ParaSystem.BShielding = ParaSystem.Bshielding_kp1;
%     end
    

    
 %% Test Shielding Field   
 
    
    B_DipCAN = MagField_DipCAN(rSol, ParaGrid, ParaSystem);
    
    rWhole = vertcat(rSol, flipud( rSol(1:end-1, :)) );
    PhiListWhole = linspace(-pi/2, 3*pi/2, size(rWhole, 1));
    [PhiGrid, ThetaGrid] = ndgrid( PhiListWhole, ParaGrid.ThetaList );
    B_DipCAN.XValWhole = B_DipCAN.XInterp(PhiGrid, ThetaGrid);
    B_DipCAN.YValWhole = B_DipCAN.YInterp(PhiGrid, ThetaGrid);
    B_DipCAN.ZValWhole = B_DipCAN.ZInterp(PhiGrid, ThetaGrid);


    if isequal([ParaSystem.BShielding.XVal, ParaSystem.BShielding.YVal, ParaSystem.BShielding.ZVal], [0, 0, 0])
        BField_km1 = B_DipCAN ;
    else
        BShield_km1 = ParaSystem.BShielding ;
        BShield_X = BShield_km1.XInterp(PhiGrid, ThetaGrid);
        BShield_Y = BShield_km1.YInterp(PhiGrid, ThetaGrid);
        BShield_Z = BShield_km1.ZInterp(PhiGrid, ThetaGrid);
        
%         BShield_X = BShield_km1.XValWhole;
%         BShield_Y = BShield_km1.YValWhole;
%         BShield_Z = BShield_km1.ZValWhole;
        
        BField_km1.XValWhole = B_DipCAN.XValWhole + BShield_X;
        BField_km1.YValWhole = B_DipCAN.YValWhole + BShield_Y;
        BField_km1.ZValWhole = B_DipCAN.ZValWhole + BShield_Z;
        BField_km1.XInterp = griddedInterpolant(PhiGrid, ThetaGrid, BField_km1.XValWhole);
        BField_km1.YInterp = griddedInterpolant(PhiGrid, ThetaGrid, BField_km1.YValWhole);
        BField_km1.ZInterp = griddedInterpolant(PhiGrid, ThetaGrid, BField_km1.ZValWhole);
    end
    
    %{
    % Visualisation
    figure;
    pcolor( B_DipCAN.ZValWhole - BField_km1.ZValWhole)
    colorbar
    %}
    

        % Defining M-Grid
%         XList = (0:1/2:ParaSystem.Nose.KSMU_Rp(1)) .';
%         ZList = (-35:1/2:35) .';
        XList = (0:1/2:40) .';
        ZList = (-35:1/2:40) .';
        [XGrid, ZGrid] = ndgrid( XList, ZList );

        BM_X_List = 0.* XGrid;
        BM_Y_List = 0.* XGrid;
        BM_Z_List = 0.* XGrid;

        CPhi = 1;
        CTheta = 1;
        for kMX=1:length(XList)
            100*(kMX/length(XList))

                for kMZ=1:length(ZList)

                    
                    M.X = XList(kMX);
                    M.Y = ZList(kMZ).*tan(ParaSystem.AlphaDipole);
                    M.Z = ZList(kMZ);

                    %{
                    dBpX = @(X) dBpM(rSol, [X(:,1), X(:,2)], [M.X, M.Y, M.Z], ParaGrid, ParaSystem, BField_km1, 'X') ;  
                    dBpY = @(X) dBpM(rSol, [X(:,1), X(:,2)], [M.X, M.Y, M.Z], ParaGrid, ParaSystem, BField_km1, 'Y') ;  
                    dBpZ = @(X) dBpM(rSol, [X(:,1), X(:,2)], [M.X, M.Y, M.Z], ParaGrid, ParaSystem, BField_km1, 'Z') ;  

                    BM_Z_List(kMX, kMY, kMZ) = integralN_mc(dBpX, [-pi/2, 3*pi/2; 0, ParaGrid.ThetaList(end)], 'timelimit', 0.1) ;
                    BM_X_List(kMX, kMY, kMZ) = integralN_mc(dBpY, [-pi/2, 3*pi/2; 0, ParaGrid.ThetaList(end)], 'timelimit', 0.1) ;
                    BM_Y_List(kMX, kMY, kMZ) = integralN_mc(dBpZ, [-pi/2, 3*pi/2; 0, ParaGrid.ThetaList(end)], 'timelimit', 0.1) ;
                    %}


                   dBp = @(X) dBpM( rSol, [X(:,1), X(:,2)], [M.X, M.Y, M.Z], ParaGrid, ParaSystem, BField_km1 ) ;  

                    PhiP = PhiGrid;
                    ThetaP = ThetaGrid;

                    dBpVal = dBp( [PhiP(:), ThetaP(:)] );

                    %{
                BM_X_List(kMX, kMY, kMZ) = sum( dBpVal(:,1) );
                BM_Y_List(kMX, kMY, kMZ) = sum( dBpVal(:,2) );
                BM_Z_List(kMX, kMY, kMZ) = sum( dBpVal(:,3) );  
                    %}              

                    % TEST: Trapezoidal Rule
                    dBpX_Grid = reshape( dBpVal(:, 1), size(PhiGrid) );
                    dBpY_Grid = reshape( dBpVal(:, 2), size(PhiGrid) );
                    dBpZ_Grid = reshape( dBpVal(:, 3), size(PhiGrid) );


                    BM_X_List(kMX, kMZ) = (1/4) .* ( dBpX_Grid(1, 1) + dBpX_Grid(end, 1) + dBpX_Grid(1, end) + dBpX_Grid(end, end) + ...
                               2 .* sum(dBpX_Grid(2:end-1, 1)) + 2 .* sum(dBpX_Grid(2:end-1, end)) + 2 .* sum(dBpX_Grid(1, 2:end-1)) + 2 .* sum(dBpX_Grid(end, 2:end-1)) + ...
                               4 .* sum( sum(dBpX_Grid(2:end-1, 2:end-1), 1)) );

                    BM_Y_List(kMX, kMZ) = (1/4) .* ( dBpY_Grid(1, 1) + dBpY_Grid(end, 1) + dBpY_Grid(1, end) + dBpY_Grid(end, end) + ...
                               2 .* sum(dBpY_Grid(2:end-1, 1)) + 2 .* sum(dBpY_Grid(2:end-1, end)) + 2 .* sum(dBpY_Grid(1, 2:end-1)) + 2 .* sum(dBpY_Grid(end, 2:end-1)) + ...
                               4 .* sum( sum(dBpY_Grid(2:end-1, 2:end-1), 1)) );

                    BM_Z_List(kMX, kMZ) = (1/4) .* ( dBpZ_Grid(1, 1) + dBpZ_Grid(end, 1) + dBpZ_Grid(1, end) + dBpZ_Grid(end, end) + ...
                               2 .* sum(dBpZ_Grid(2:end-1, 1)) + 2 .* sum(dBpZ_Grid(2:end-1, end)) + 2 .* sum(dBpZ_Grid(1, 2:end-1)) + 2 .* sum(dBpZ_Grid(end, 2:end-1)) + ...
                               4 .* sum( sum(dBpZ_Grid(2:end-1, 2:end-1), 1) ) );
            %{               
            dBp = @(X) dBpM(rSol, [X(:,1), X(:,2)], [M.X, M.Y, M.Z], ParaGrid, ParaSystem, BField_km1, 'Z') ;
            integralN_mc(dBp, [-pi/2, 3*pi/2; 0, ParaGrid.ThetaList(end)], 'timelimit', 1) 
            %}
                               
            
        % TEST: integral2
%             BxInterp = griddedInterpolant(PhiGrid, ThetaGrid, dBpX_Grid, 'spline');
%                 BM_X_List(kMX, kMZ) = quad2d( @(phi,theta) BxInterp(phi,theta), -pi/2, 3*pi/2, 0, ParaGrid.ThetaList(end) );
%             ByInterp = griddedInterpolant(PhiGrid, ThetaGrid, dBpY_Grid, 'spline');
%                 BM_Y_List(kMX, kMZ) = quad2d( @(phi,theta) ByInterp(phi,theta), -pi/2, 3*pi/2, 0, ParaGrid.ThetaList(end) );
%             BzInterp = griddedInterpolant(PhiGrid, ThetaGrid, dBpZ_Grid, 'spline');
%                 BM_Z_List(kMX, kMZ) = quad2d( @(phi,theta) BzInterp(phi,theta), -pi/2, 3*pi/2, 0, ParaGrid.ThetaList(end) );
        % TEST: integral2

                                %{
                               % Visualisation
                               figure;
                               pcolor(dBpZ_Grid)
                               colorbar
                                %}
                end
        end


    
    % Visualisation: Extracting Fields
    [XGrid_Val, YGrid_Val, ZGrid_Val] = ndgrid(XList, ZList.*tan(ParaSystem.AlphaDipole), ZList);
    Grid.X = reshape(XGrid_Val(:,1,:), [length(XList), length(ZList)]);
    Grid.Y = reshape(YGrid_Val(:,1,:), [length(XList), length(ZList)]);
    Grid.Z = reshape(ZGrid_Val(:,1,:), [length(XList), length(ZList)]);
    BDipole = DipoleField(Grid, ParaSystem);
        
    BCAN_NMM.XVal = BCAN_Rotated.XInterp(XGrid_Val, YGrid_Val, ZGrid_Val);
    BCAN_NMM.X = reshape( BCAN_NMM.XVal(:,1,:), [length(XList), length(ZList)] );
    BCAN_NMM.YVal = BCAN_Rotated.YInterp(XGrid_Val, YGrid_Val, ZGrid_Val);
    BCAN_NMM.Y = reshape( BCAN_NMM.YVal(:,1,:), [length(XList), length(ZList)] );
    BCAN_NMM.ZVal = BCAN_Rotated.ZInterp(XGrid_Val, YGrid_Val, ZGrid_Val);
    BCAN_NMM.Z = reshape( BCAN_NMM.ZVal(:,1,:), [length(XList), length(ZList)] );
    
%     If Phi > zero
    %{
    [XGrid_Val, YGrid_Val, ZGrid_Val] = ndgrid(XList, ZList.*tan(ParaSystem.AlphaDipole), ZList);
    Grid.X = reshape(XGrid_Val(:,1,:), [length(XList), length(ZList)]);
    Grid.Z = reshape(ZGrid_Val(:,1,:), [length(XList), length(ZList)]);
    Grid.Y = Grid.Z .* tan(ParaSystem.AlphaDipole);
    BDipole = DipoleField(Grid, ParaSystem);

    
    BCAN_NMM.XVal = BCAN_Rotated.XInterp(XGrid_Val, ZGrid_Val.*tan(ParaSystem.AlphaDipole), ZGrid_Val);
    BCAN_NMM.X = reshape( BCAN_NMM.XVal(:,1,:), [length(XList), length(ZList)] );
    BCAN_NMM.YVal = BCAN_Rotated.YInterp(XGrid_Val, ZGrid_Val.*tan(ParaSystem.AlphaDipole), ZGrid_Val);
    BCAN_NMM.Y = reshape( BCAN_NMM.YVal(:,1,:), [length(XList), length(ZList)] );
    BCAN_NMM.ZVal = BCAN_Rotated.ZInterp(XGrid_Val, ZGrid_Val.*tan(ParaSystem.AlphaDipole), ZGrid_Val);
    BCAN_NMM.Z = reshape( BCAN_NMM.ZVal(:,1,:), [length(XList), length(ZList)] );
    
    figure;
    hold on
        S = [1, 1, 1];
        B_NMM = @(S) Btot.Z(S) .* cos(ParaSystem.AlphaDipole) +  Btot.Y(S) .* sin(ParaSystem.AlphaDipole);
        Z_Alpha = Z_NMM .* cos(ParaSystem.AlphaDipole);
        SL1 = streamline(stream2( (X_NMM).', (Z_Alpha).', (Btot.X(S)).', (B_NMM(S)).', XStartList, 0.*ZStartList, [0.1 1000]));
        SL2 = streamline(stream2( (X_NMM).', (Z_Alpha).', -(Btot.X(S)).', -(B_NMM(S)).', XStartList, 0.*ZStartList, [0.1 1000]));
        axis equal

    %}
        
    BShield.X = BM_X_List;  
    BShield.Y = BM_Y_List;  
    BShield.Z = BM_Z_List;  

    % Visualisation: Defining NMM Fields
    X_NMM = reshape(XGrid_Val(:,1,:), [length(XList), length(ZList)]);
    Y_NMM = reshape(YGrid_Val(:,1,:), [length(XList), length(ZList)]);
    Z_NMM = reshape(ZGrid_Val(:,1,:), [length(XList), length(ZList)]);
    
    Btot.X = @(S) S(1) .* BDipole.X + S(2) .* BCAN_NMM.X + S(3) .* BShield.X;
    Btot.Y = @(S) S(1) .* BDipole.Y + S(2) .* BCAN_NMM.Y + S(3) .* BShield.Y;
    Btot.Z = @(S) S(1) .* BDipole.Z + S(2) .* BCAN_NMM.Z + S(3) .* BShield.Z;
    Btot.Norm = @(S) sqrt(Btot.X(S).^2 + Btot.Y(S).^2 + Btot.Z(S).^2 );

    
% Plotting: Magnetodisk Field
figure;
hold on
    S = [1, 1, 0];
%     XStartList = linspace(0, 35, 30);
%     ZStartList = XStartList .* tan(-27*pi/180);
    XStartList = linspace(0, 35, 20);
%     ZStartList = XStartList .* tan( -acos(ParaSystem.NoseDirection(1)) );
    ZStartList = XStartList .* tan(-27*pi/180)  ;
    
    SL1 = streamline(stream2( (X_NMM).', (Z_NMM).', (Btot.X(S)).', (Btot.Z(S)).', XStartList, ZStartList, [0.1 2000]));
    SL2 = streamline(stream2( (X_NMM).', (Z_NMM).', -(Btot.X(S)).', -(Btot.Z(S)).', XStartList, ZStartList, [0.1 2000]));

    
    ColorSL = cbrewer2('Purples', 20);
    Opacity = 0.5;
    set([SL1, SL2], 'LineWidth', 3, 'Color', [ColorSL(17,:), Opacity])
    
%     Equator_Rot = line( [XStartList(1), 27], [XStartList(1) .* tan(-27*pi/180), 27 .* tan(-27*pi/180)] );
    Equator_Rot = line( [XStartList(1), 23], [XStartList(1) .* tan(-27*pi/180), 23 .* tan(-27*pi/180) ] );
    
    ColorEquRot = cbrewer2('Blues', 20);
    Opacity = 0.2;
    set(Equator_Rot, 'LineStyle', '-', 'LineWidth', 7, 'Color', [ColorEquRot(19,:), Opacity] );

    % Axes
    axis equal
%     axis([-3 30 -32 10]);
    axis([-3 35 -20 20]);

    axes = gca;
    set(axes, 'FontName'   , 'Helvetica',...
        'XAxisLocation', 'top',...
        'YAxisLocation', 'right',...
        'FontSize', 16,...
        'LineWidth'   , 1 ...
        );
    Xlabel = xlabel('X_0 (R_p)');
    Ylabel = ylabel('Z_0 (R_p)');
    grid on
    ColourBackground = [0, 0, 0, 0.3];
    set(gca,'Color',ColourBackground)
    set(gcf,'color','w');
    Xaxis = line([-5 40], [0 0])   ;
    Zaxis = line([0 0], [- 35 35])   ;
    set([Xaxis, Zaxis], 'LineWidth', 2, 'LineStyle', '--', 'Color', [0 0 0 0.3])   
                
    % Plotting: Dipole Field
    S = [1, 0, 0];
    XStartList = linspace(0, 35, 20);
%     ZStartList = XStartList .* tan(-27*pi/180);
    ZStartList = XStartList .* tan(-27*pi/180).*0;
    SL3 = streamline(stream2( (X_NMM).', (Z_NMM).', (Btot.X(S)).', (Btot.Z(S)).', XStartList, ZStartList, [0.1 2000]));
    SL4 = streamline(stream2( (X_NMM).', (Z_NMM).', -(Btot.X(S)).', -(Btot.Z(S)).', XStartList, ZStartList, [0.1 2000]));
    
    ColorSL = cbrewer2('Greens', 20);
    Opacity2 = 0.3;
    set([SL3, SL4], 'LineWidth', 2, 'Color', [ColorSL(19,:), Opacity2])
    
    set (gca,'Xdir','reverse')

    
export_fig FieldLines_Magnetodisk_Summer25.pdf -opengl -q101
    
    
    
    % Plotting: Total Field
figure;   
hold on
    S = [1, 1, 1];
%     XStartList = linspace(0, 35, 30);
%     ZStartList = XStartList .* tan(-27*pi/180);
    XStartList = linspace(0, 22, 30);
%     ZStartList = XStartList .* tan( -acos(ParaSystem.NoseDirection(1)) );
    ZStartList = XStartList .* tan(-27*pi/180);
    
%     scatter3(XStartList, ZStartList, 0.*ZStartList, 'filled')
    SL1 = streamline(stream2( (X_NMM).', (Z_NMM).', (Btot.X(S)).', (Btot.Z(S)).', XStartList, ZStartList, [0.1 2000]));
    SL2 = streamline(stream2( (X_NMM).', (Z_NMM).', -(Btot.X(S)).', -(Btot.Z(S)).', XStartList, ZStartList, [0.1 2000]));
    
    ColorSL = cbrewer2('PuRd', 20);
    Opacity = 0.3;
    set([SL1, SL2], 'LineWidth', 3, 'Color', [ColorSL(19,:), Opacity]);
    
%     Equator_Rot = line( [XStartList(1), XStartList(end)], [XStartList(1) .* tan(-27*pi/180), XStartList(end) .* tan(-27*pi/180)] );
%     ColorEquRot = cbrewer2('Blues', 20);
%     Opacity = 0.5;
%     set(Equator_Rot, 'LineStyle', ':', 'LineWidth', 3, 'Color', [ColorEquRot(17,:), Opacity] );

    % Axes
    axis equal
    axis([-3 30 -32 10]);
    axes = gca;
    set(axes, 'FontName'   , 'Helvetica',...
        'XAxisLocation', 'top',...
        'YAxisLocation', 'right',...
        'FontSize', 16,...
        'LineWidth'   , 1 ...
        );
    Xlabel = xlabel('X_0 (R_p)');
    Ylabel = ylabel('Z_0 (R_p)');
    grid on
    ColourBackground = [0, 0, 0, 0.3];
    set(gca,'Color',ColourBackground)
    set(gcf,'color','w');
    Xaxis = line([-5 40], [0 0])   ;
    Zaxis = line([0 0], [- 35 35])   ;
    set([Xaxis, Zaxis], 'LineWidth', 2, 'LineStyle', '--', 'Color', [0 0 0 0.3])   

    % Find Magnetic Equator    
    EquationSL_1 = stream2( X_NMM.', Z_NMM.', -Btot.X(S).',-Btot.Z(S).', XStartList, ZStartList);
    RhoReversal.XUp = zeros(length(EquationSL_1), 1);
    RhoReversal.ZUp = zeros(length(EquationSL_1), 1);    
    for kLine = 1:length(EquationSL_1)
        Array = EquationSL_1{1,kLine};
        if not(isequal(Array, []))
            X = Array(:, 1);
            X = X(not(isnan(X)));
            Z = Array(:, 2);
            Z = Z(not(isnan(Z)));
        
            Position = [X, Z];
            eEqu = [cos(-27*pi/180); sin(-27*pi/180)];
            rho = Position(:,1) .* eEqu(1) + Position(:,2) .* eEqu(2);
%             rhoMax = max(rho);
            rho = sqrt(X.^2+Z.^2);
            rhoMax = max(rho);
        
            XMax= X(rho == rhoMax);
            ZCorr = Z(X == XMax);
            RhoReversal.XUp(kLine) = XMax;
            RhoReversal.ZUp(kLine) = ZCorr;
        end
    end
    
    Equator_Mag = plot(RhoReversal.XUp, RhoReversal.ZUp);        
    ColorEquMag = cbrewer2('Blues', 20);
    Opacity = 0.2;
    set(Equator_Mag, 'LineStyle', '-', 'LineWidth', 7, 'Color', [ColorEquMag(19,:), Opacity] );
    
     set (gca,'Xdir','reverse')
     
     Nose = scatter3(ParaSystem.Nose.KSMU_Rp(1), ParaSystem.Nose.KSMU_Rp(3), 0, 'filled');
     ColorNose = cbrewer2('Reds', 20);
     Opacity = 0.7;
     set(Nose, 'MarkerFaceColor', [ColorNose(18,:)], 'SizeData', 150, 'MarkerFaceAlpha', Opacity);
     
     
 export_fig FieldLines_TotalField_Summer25_Arridge.pdf -opengl -q101
    
 
% TEST: Arridge Bowl
x = (0:26);
z = 0.*x ;
r = sqrt(x.^2+z.^2);

Rh = 29;
ThetaSun = 27.*pi/180 ;

AngleRot = ThetaSun;
 
Coord = [x.', z.'];
RotCoord_X = cos(-AngleRot) .* Coord(:,1) - sin(-AngleRot) .* Coord(:,2);
RotCoord_Z = sin(-AngleRot) .* Coord(:,1) + cos(-AngleRot) .* Coord(:,2);

Cond = r>=0;
Offset = ( r - Rh.*tanh(r./Rh) ) .* tan(ThetaSun) ;

z_cs = RotCoord_Z;
x_cs = RotCoord_X;
z_cs(Cond) = RotCoord_Z(Cond) - cos(AngleRot).*Offset(Cond) .';
x_cs(Cond) = RotCoord_X(Cond) - sin(AngleRot).*Offset(Cond) .';

hold on
    PlotArridge = plot(x_cs, z_cs);
    ColorCS = cbrewer2('YlOrBr', 20);
    Opacity = 0.3;
    set(PlotArridge, 'LineWidth', 3, 'Color', [ColorCS(15,:), Opacity], 'LineStyle', '-.');
    
%     scatter(Coord(:,1), Coord(:,2))
%     scatter(RotCoord_X, RotCoord_Z)
    axis equal
 
% TEST: Arridge Bowl
 

%% TEST: Add MP surface
P_CartNose2KSM = ParaSystem.P_CartNose2KSM;

PhiList = [-pi/2; pi/2];
ThetaList = linspace(ParaGrid.ThetaList(1), ParaGrid.ThetaList(end), length(ParaGrid.ThetaList)*2);
[PhiGrid, ThetaGrid] = ndgrid(PhiList, ThetaList);
rNMM = rInterp(PhiGrid, ThetaGrid);
[Y, Z, X] = sph2cart(PhiGrid, pi/2-ThetaGrid, rNMM.*ParaSystem.r0/ParaSystem.Rp);
X_KSM = P_CartNose2KSM(1,1) .* X + P_CartNose2KSM(1,2) .* Y + P_CartNose2KSM(1,3) .* Z;
Y_KSM = P_CartNose2KSM(2,1) .* X + P_CartNose2KSM(2,2) .* Y + P_CartNose2KSM(2,3) .* Z;
Z_KSM = P_CartNose2KSM(3,1) .* X + P_CartNose2KSM(3,2) .* Y + P_CartNose2KSM(3,3) .* Z;

% figure;
hold on
    MP_South = plot(X_KSM(1,:), Z_KSM(1,:));
    MP_North = plot(X_KSM(2,:), Z_KSM(2,:));
    ColorPlot = cbrewer2('YlGn', 20);
    Opacity = 0.7;
    set([MP_South, MP_North], 'linewidth', 7, 'linestyle', '-.', 'Color', [ColorPlot(18,:), Opacity])
    axis equal
    
    axis([-3 35 -28 10])

%% TEST: Add MP surface





    
S = [1; 1; 1]
 figure
    hold on
        PlotQuiver2 = quiver(X_NMM, Z_NMM, Btot.X(S) ./Btot.Norm(S), Btot.Z(S) ./Btot.Norm(S), 0);
        axis equal
        
    XStartList = linspace(0, 28, 30);
    ZStartList = XStartList .* tan(-27*pi/180);
    Cond = sqrt(X_NMM.^2 + Z_NMM.^2) > 5;
    SL1 = streamline(stream2( (X_NMM).', (Z_NMM).', (Btot.X).', (Btot.Z).', XStartList, ZStartList, [0.1 1000]));
    SL2 = streamline(stream2( (X_NMM).', (Z_NMM).', -(Btot.X).', -(Btot.Z).', XStartList, ZStartList, [0.1 1000]));
        
        
    
        
    % Find Magnetic Equator    
    EquationSL_1 = stream2( X_NMM.', Z_NMM.', -Btot.X.',-Btot.Z.', XStartList, ZStartList);
    RhoReversal.XUp = zeros(length(EquationSL_1), 1);
    RhoReversal.ZUp = zeros(length(EquationSL_1), 1);    
    for kLine = 1:length(EquationSL_1)
        Array = EquationSL_1{1,kLine};
        if not(isequal(Array, []))
            X = Array(:, 1);
            X = X(not(isnan(X)));
            Z = Array(:, 2);
            Z = Z(not(isnan(Z)));
        
            Position = [X, Z];
            eEqu = [cos(-27*pi/180); sin(-27*pi/180)];
            rho = Position(:,1) .* eEqu(1) + Position(:,2) .* eEqu(2);
            rhoMax = max(rho);
        
            XMax= X(rho == rhoMax);
            ZCorr = Z(X == XMax);
            RhoReversal.XUp(kLine) = XMax;
            RhoReversal.ZUp(kLine) = ZCorr;
        end
    end
    
    Plot_RhoReversal = scatter(RhoReversal.XUp, RhoReversal.ZUp);        
        
        

        
        

    
    
    
    
    
    
    
    
    
    
    
    
    Method = 'linear';
    
    BShield.XInterp = griddedInterpolant(XGrid, YGrid, ZGrid, BM_X_List, Method);
    BShield.YInterp = griddedInterpolant(XGrid, YGrid, ZGrid, BM_Y_List, Method);
    BShield.ZInterp = griddedInterpolant(XGrid, YGrid, ZGrid, BM_Z_List, Method);

    XVal = (0:1:40) .';
%     YVal = (0:1:35).';
    YVal = 0;
    ZVal = (-35:1:35).';
    [XGridVal, YGridVal, ZGridVal] = ndgrid(XVal, YVal, ZVal);
    
    % Visualisation
    Bshield_X = BShield.XInterp(XGridVal, YGridVal, ZGridVal);
    Bshield_Y = BShield.YInterp(XGridVal, YGridVal, ZGridVal);
    Bshield_Z = BShield.ZInterp(XGridVal, YGridVal, ZGridVal);
    
    BCAN_X = BCAN_Rotated.XInterp(XGridVal, YGridVal, ZGridVal);
    BCAN_Y = BCAN_Rotated.YInterp(XGridVal, YGridVal, ZGridVal);
    BCAN_Z = BCAN_Rotated.ZInterp(XGridVal, YGridVal, ZGridVal);
    
    Grid.X = XGridVal;
    Grid.Y = YGridVal;
    Grid.Z = ZGridVal;
    BDip = DipoleField(Grid, ParaSystem);
    BDipole_X =  BDip.X;
    BDipole_Y =  BDip.Y;
    BDipole_Z =  BDip.Z;

    S.Dip = 1;
    S.CAN = 1;
    S.Shield = 1;
    Btot.X = S.Dip .* BDipole_X + S.CAN .* BCAN_X + S.Shield .* Bshield_X ;
    Btot.Y = S.Dip .* BDipole_Y + S.CAN .* BCAN_Y + S.Shield .* Bshield_Y ;
    Btot.Z = S.Dip .* BDipole_Z + S.CAN .* BCAN_Z + S.Shield .* Bshield_Z ;
    Btot.Norm = sqrt( Btot.X.^2 + Btot.Y.^2 + Btot.Z.^2 );

    Cond  = sqrt(XGridVal.^2 + YGridVal.^2 + ZGridVal.^2) > 5;
    figure;
    hold on
        PlotQuiver2 = quiver3(XGridVal(Cond), YGridVal(Cond), ZGridVal(Cond), Btot.X(Cond) ./ Btot.Norm(Cond) , Btot.Y(Cond) ./ Btot.Norm(Cond), Btot.Z(Cond) ./ Btot.Norm(Cond), 0);
%         PlotQuiver3 = quiver3(XGridVal(Cond), YGridVal(Cond), ZGridVal(Cond), BDipole_X(Cond)+BCAN_X(Cond), BDipole_Y(Cond)+BCAN_Y(Cond), BDipole_Z(Cond)+BCAN_Z(Cond), 0);
        set( [PlotQuiver2], 'AutoScaleFactor', 1, 'LineWidth', 1, 'MaxHeadSize', 1)
    line( [0 40], [0 0], [0 -40*tan(27*pi/180)] )
    axis equal
    axis([0 35 -1/2 1/2 -20 5])
    
    
    
    
% --- Testing Plotting Streamlines
    XVal = (0:1:35) .';
    YVal = (0:1:35).';
    ZVal = (-35:1:35).';
    [XGridVal_Mesh, YGridVal_Mesh, ZGridVal_Mesh] = meshgrid(XVal, YVal, ZVal);
    
    NLines = 40;
    XStartList = linspace(7, 30, NLines );
    ZStartList = XStartList .* tan(-27*pi/180);

    Bshield_X = BShield.XInterp(XGridVal_Mesh, YGridVal_Mesh, ZGridVal_Mesh);
    Bshield_Y = BShield.YInterp(XGridVal_Mesh, YGridVal_Mesh, ZGridVal_Mesh);
    Bshield_Z = BShield.ZInterp(XGridVal_Mesh, YGridVal_Mesh, ZGridVal_Mesh);
    
    BCAN_X = BCAN_Rotated.XInterp(XGridVal_Mesh, YGridVal_Mesh, ZGridVal_Mesh);
    BCAN_Y = BCAN_Rotated.YInterp(XGridVal_Mesh, YGridVal_Mesh, ZGridVal_Mesh);
    BCAN_Z = BCAN_Rotated.ZInterp(XGridVal_Mesh, YGridVal_Mesh, ZGridVal_Mesh);
    
    Grid.X = XGridVal_Mesh;
    Grid.Y = YGridVal_Mesh;
    Grid.Z = ZGridVal_Mesh;
    BDip = DipoleField(Grid, ParaSystem);
    BDipole_X =  BDip.X;
    BDipole_Y =  BDip.Y;
    BDipole_Z =  BDip.Z;
    
    S.Dip = 1;
    S.CAN = 1;
    S.Shield = 1;
    
    Btot.X = S.Dip .* BDipole_X + S.CAN .* BCAN_X + S.Shield .* Bshield_X ;
    Btot.Y = S.Dip .* BDipole_Y + S.CAN .* BCAN_Y + S.Shield .* Bshield_Y ;
    Btot.Z = S.Dip .* BDipole_Z + S.CAN .* BCAN_Z + S.Shield .* Bshield_Z ;

    Btot_NMM_X = reshape(Btot.X(YGridVal_Mesh==0), [length(XVal), length(ZVal)]);
    Btot_NMM_Y = reshape(Btot.Y(YGridVal_Mesh==0), [length(XVal), length(ZVal)]);
    Btot_NMM_Z = reshape(Btot.Z(YGridVal_Mesh==0), [length(XVal), length(ZVal)]);
    
    X_NMM = reshape(XGridVal_Mesh(YGridVal_Mesh==0), [length(XVal), length(ZVal)]);
    Z_NMM = reshape(ZGridVal_Mesh(YGridVal_Mesh==0), [length(XVal), length(ZVal)]);

    % Plotting Field Lines
    C = 2;
    figure;
    hold on
    Plot_StreamlinesUp = streamline( X_NMM(1:C:end, 1:C:end).', Z_NMM(1:C:end, 1:C:end).', Btot_NMM_X(1:C:end, 1:C:end).', Btot_NMM_Z(1:C:end, 1:C:end).', XStartList(1:C:end, 1:C:end), ZStartList(1:C:end, 1:C:end), [0.01, 20000]);
    Plot_StreamlinesDown = streamline( X_NMM(1:C:end, 1:C:end).', Z_NMM(1:C:end, 1:C:end).', -Btot_NMM_X(1:C:end, 1:C:end).', -Btot_NMM_Z(1:C:end, 1:C:end).', XStartList(1:C:end, 1:C:end), ZStartList(1:C:end, 1:C:end), [0.01, 20000]);
    RotEquator = line( [0 40], [0 -40*tan(27*pi/180)] ) ;
    axis equal
    
    % Finding Displaced Equator
    EquationSL_1 = stream2( X_NMM(1:C:end, 1:C:end).', Z_NMM(1:C:end, 1:C:end).', Btot_NMM_X(1:C:end, 1:C:end).', Btot_NMM_Z(1:C:end, 1:C:end).', XStartList(1:C:end, 1:C:end), ZStartList(1:C:end, 1:C:end));
    RhoReversal.XUp = zeros(length(EquationSL_1), 1);
    RhoReversal.ZUp = zeros(length(EquationSL_1), 1);    
    for kLine = 1:length(EquationSL_1)
        Array = EquationSL_1{1,kLine};
        if not(isequal(Array, []))
            X = Array(:, 1);
            X = X(not(isnan(X)));
            Z = Array(:, 2);
            Z = Z(not(isnan(Z)));
            XMax = max(X);
            ZCorr = Z(X == XMax);
            RhoReversal.XUp(kLine) = XMax;
            RhoReversal.ZUp(kLine) = ZCorr;
        end
    end
    
%     Plot_RhoReversal = scatter(RhoReversal.XUp, RhoReversal.ZUp);
    
    % Options
    ColorSL = cbrewer2('Greys', 20);
    ColorReversal = cbrewer2('Blues', 20);
    Opacity = 0.5;
    set([Plot_StreamlinesUp, Plot_StreamlinesDown], 'LineWidth', 2, 'Color', horzcat( ColorSL(17, :), Opacity) ) ;
    % set(RotEquator, 'LineWidth', 2, 'LineStyle', '--')
%     set(Plot_RhoReversal, 'Marker', 'x', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', ColorReversal(17, :), 'MarkerFaceAlpha', 0.5, 'SizeData', 50);
     set(Plot_RhoReversal, 'Marker', '+', 'LineWidth', 2, 'MarkerFaceColor', ColorReversal(19, :), 'MarkerFaceAlpha', 0.3, 'SizeData', 50)
    
    % Axes
    axes = gca;
    set(axes, 'FontName'   , 'Helvetica',...
        'XAxisLocation', 'bottom',...
        'YAxisLocation', 'left',...
        'FontSize', 16,...
        'LineWidth'   , 1 ...
        );
    Xlabel = xlabel('X (R_p)');
    Ylabel = ylabel('Z (R_p)');
    grid on
    ColourBackground = [0, 0, 0, 0.3];
    set(gca,'Color',ColourBackground)
    set(gcf,'color','w');
    axis([-3 37 -27 10]);
    Xaxis = line([-5 40], [0 0])   ;
    Zaxis = line([0 0], [- 35 35])   ;
    set([Xaxis, Zaxis], 'LineWidth', 2, 'LineStyle', '--', 'Color', [0 0 0 0.3])   
        
    scatter( X_KSM( abs(ParaGrid.PhiGrid)==pi/2 ), Z_KSM( abs(ParaGrid.PhiGrid)==pi/2 ), 'filled' )


    
    export_fig FieldLinesNMM.pdf -opengl -q101
    

    
    
        
    
    XNose = P_KSM2CartNose(1,1) .* XGridVal +  P_KSM2CartNose(1,2) .* YGridVal + P_KSM2CartNose(1,3) .* ZGridVal ;
    YNose = P_KSM2CartNose(2,1) .* XGridVal +  P_KSM2CartNose(2,2) .* YGridVal + P_KSM2CartNose(2,3) .* ZGridVal ;
    ZNose = P_KSM2CartNose(3,1) .* XGridVal +  P_KSM2CartNose(3,2) .* YGridVal + P_KSM2CartNose(3,3) .* ZGridVal ;
    
    [Phi, Elevation, r] = cart2sph( YNose, ZNose, XNose );
    Theta = pi/2 - Elevation;
    BCANX = BDipCAN.XInterp(Phi, Theta);
    BCANY = BDipCAN.YInterp(Phi, Theta);
    BCANZ = BDipCAN.ZInterp(Phi, Theta);
    
    figure;
    hold on
    PlotQuiver1 = quiver3( XGridVal, YGridVal, ZGridVal, BX.Val + BCANX, BY.Val + BCANY, BZ.Val + BCANZ );
    PlotQuiver2 = quiver3( XGridVal, YGridVal, ZGridVal, BCANX, BCANY, BCANZ );
    line( [0 40], [0 0], [0 -40*tan(27*pi/180)] )
    set( [PlotQuiver1, PlotQuiver2], 'AutoScaleFactor', 1, 'LineWidth', 1, 'MaxHeadSize', 1)
    axis equal
    axis( [0 35 0 35 -35 0] )
    
    
    XVal = (-0:2:35).';
    YVal = (-0:2:35).';
    ZVal = (-35:2:0).';
    [XGridVal, YGridVal, ZGridVal] = meshgrid(XVal, YVal, ZVal);
    
    XStart = XVal;
    YStart = 0.*YVal;
    ZStart = -XVal.*tan(27*pi/180);
    U = BCANX + BX.Val;
    V = BCANY + BY.Val;
    W = BCANZ + BZ.Val;
    figure;
    hold on
    streamline( XGridVal, YGridVal, ZGridVal, U, V, W, XStart, YStart, ZStart )
    
    
    
%% -------------------------------------------------------------------------------
% Visualisation of Optimised Surface

% ----  Computation of DeltaP/MeanP Metric
    ResidualGrid_Sol = log10( abs( PB_Grid(rSol, ParaGrid, ParaSystem) ) );

% ----  Building Optimised Surface: FIRST Half
    Elevation = pi/2 - ParaGrid.ThetaGrid;
    [XNoseMB, YNoseMB, ZNoseMB] = sph2cart(ParaGrid.PhiGrid, Elevation, rSol .* ParaSystem.r0 ./ ParaSystem.Rp);
    YNose = XNoseMB;
    ZNose = YNoseMB;
    XNose = ZNoseMB;
    P_CartNose2KSM = ParaSystem.P_CartNose2KSM;
    X_KSM = P_CartNose2KSM(1,1) .* XNose + P_CartNose2KSM(1,2) .* YNose + P_CartNose2KSM(1,3) .* ZNose;
    Y_KSM = P_CartNose2KSM(2,1) .* XNose + P_CartNose2KSM(2,2) .* YNose + P_CartNose2KSM(2,3) .* ZNose;
    Z_KSM = P_CartNose2KSM(3,1) .* XNose + P_CartNose2KSM(3,2) .* YNose + P_CartNose2KSM(3,3) .* ZNose;
    
    figure;
    hold on
    Surface = surf(X_KSM, Y_KSM, Z_KSM, X_KSM);
    XValue = 16;
    Contours = plot3(X_KSM, Y_KSM, Z_KSM);

% ----  Building Optimised Surface: SECOND Half
    YNose = -YNose;           
    X_KSM = P_CartNose2KSM(1,1) .* XNose + P_CartNose2KSM(1,2) .* YNose + P_CartNose2KSM(1,3) .* ZNose;
    Y_KSM = P_CartNose2KSM(2,1) .* XNose + P_CartNose2KSM(2,2) .* YNose + P_CartNose2KSM(2,3) .* ZNose;
    Z_KSM = P_CartNose2KSM(3,1) .* XNose + P_CartNose2KSM(3,2) .* YNose + P_CartNose2KSM(3,3) .* ZNose;
    Surface2 = surf(X_KSM, Y_KSM, Z_KSM, X_KSM);
    Contours2 = plot3(X_KSM, Y_KSM, Z_KSM);
    
% ----  Plot Options
    set([Surface, Surface2], 'EdgeAlpha', 0, 'FaceAlpha', 0.6, 'FaceColor', 'texturemap', 'Cdata', ResidualGrid_Sol);
    set([Contours, Contours2], 'LineWidth', 1, 'Color', [0 0 0 0.6])

% ----  Axes / cb Options
    cb = colorbar( 'XTickLabel',{'0.1','1','10','100'}, 'XTick', [-3, -2, -1, 0] );
    axis equal
    caxis([-3 0])
    axis([-30 30 -30 30 -30 30])
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    rNMM1 = rSol( (ParaGrid.PhiGrid) == pi/2 );
    rNMM2 = rSol( (ParaGrid.PhiGrid) == -pi/2 );
    [X1,Y1] = pol2cart(ParaGrid.ThetaList, rNMM1 .* ParaSystem.r0 ./ ParaSystem.Rp);
    [X2,Y2] = pol2cart(ParaGrid.ThetaList, rNMM2 .* ParaSystem.r0 ./ ParaSystem.Rp);

    X1(end-4) = (1/2) .* (X1(end-3) + X1(end-5));
    Y1(end-4) = (1/2) .* (Y1(end-3) + Y1(end-5));
    
    Method = 'spline';
    Interp.r1 = griddedInterpolant( ParaGrid.ThetaList , rNMM1 .* ParaSystem.r0 ./ ParaSystem.Rp, Method );
    Interp.r2 = griddedInterpolant( ParaGrid.ThetaList(1:2:end), rNMM2(1:2:end) .* ParaSystem.r0 ./ ParaSystem.Rp, Method );
    ThetaVal = (ParaGrid.ThetaList(1): 3 *pi/180:ParaGrid.ThetaList(end));
    r1Val = Interp.r1(ThetaVal);
    r2Val = Interp.r2(ThetaVal);
    [XVal1,YVal1] = pol2cart(ThetaVal, r1Val);
    [XVal2,YVal2] = pol2cart(ThetaVal, r2Val);
    
    figure;
    hold on
        Plot1 = scatter( Y1 , X1);
        Plot = scatter( Y2 , X2);
        set(Plot, 'LineStyle', '--', 'LineWidth', 3);
    axis equal
    
%% Extracting NMM profiles

ThetaList = (0:1:90) .* pi/180;
PhiList = [-pi/2; 0; pi/2];
[PhiGrid, ThetaGrid] = ndgrid(PhiList, ThetaList);
rNMM = rInterp(PhiGrid, ThetaGrid);
[XNose, YNose, ZNose] = sph2cart(PhiGrid, pi/2- ThetaGrid, rNMM .* ParaSystem.r0 ./ ParaSystem.Rp);
figure;
YNose_North = YNose(3, :);
XNose_Equator = XNose(2, :);
YNose_South = YNose(1, :);
ZNose_North = ZNose(3, :);
ZNose_Equator = ZNose(2, :);
ZNose_South = ZNose(1, :);
hold on

    % Profiles
    YNose_North_2 = smoothdata(YNose_North(1:end), 2, 'sgolay');
    PlotNorth = plot(ZNose_North(:), YNose_North_2(:));
    YNose_South_2 = smoothdata(YNose_South(1:end), 3, 'sgolay');
    PlotSouth = plot( ZNose_South(:), YNose_South_2(:));
    PlotEquator_1 = plot( ZNose_Equator(:), XNose_Equator(:));
    PlotEquator_2 = plot( ZNose_Equator(:), -XNose_Equator(:));
    ColorNMM = cbrewer2('BuGn', 20);
    set([PlotNorth, PlotSouth], 'LineWidth', 7, 'Color', horzcat(ColorNMM(16,:), 0.7));
    ColorEquator = cbrewer2('YlOrBr', 10);
    set([PlotEquator_1, PlotEquator_2], 'LineStyle', '-.', 'LineWidth', 7, 'Color', horzcat(ColorEquator(9,:), 0.3))
    axis equal
    grid on

    % Disk
    a = ParaSystem.CAN_DiskParameters(2);
    b = ParaSystem.CAN_DiskParameters(3);
    D = ParaSystem.CAN_DiskParameters(4);
    V11 = [ -b; D ];
    V12 = [  b; D ];
    V21 = [ -b; -D];
    V22 = [  b; -D];
    AngleDisk = -acos(ParaSystem.NoseDirection(1)) + 27*pi/180;
    SunRotation = [ cos(AngleDisk), sin(AngleDisk); -sin(AngleDisk), cos(AngleDisk)];
    V11_Rot = SunRotation * V11;
    V12_Rot = SunRotation * V12;
    V21_Rot = SunRotation * V21;
    V22_Rot = SunRotation * V22;
    Disk1 = patch('Vertices', [V11_Rot(1), V12_Rot(1), V22_Rot(1), V21_Rot(1); V11_Rot(2), V12_Rot(2), V22_Rot(2), V21_Rot(2)].','Faces',[1 2 3 4] );
    set(Disk1, 'FaceColor', [0 0 0], 'FaceAlpha', 0.2)
    V11 = [ -a; D ];
    V12 = [  a; D ];
    V21 = [ -a; -D];
    V22 = [  a; -D];
    V11_Rot = SunRotation * V11;
    V12_Rot = SunRotation * V12;
    V21_Rot = SunRotation * V21;
    V22_Rot = SunRotation * V22;
    Disk2 = patch('Vertices', [V11_Rot(1), V12_Rot(1), V22_Rot(1), V21_Rot(1); V11_Rot(2), V12_Rot(2), V22_Rot(2), V21_Rot(2)].','Faces',[1 2 3 4] );
    ColorDisk = cbrewer2('PuBu', 10);
    set([Disk1, Disk2], 'EdgeColor', 'none', 'FaceColor', ColorDisk(9,:), 'FaceAlpha', 0.2)
    
    
    % Cusps
    CNorth = [10.12; 17.98];
    CSouth = [4.625; -26.23];
    LineCNorth = line( [0 CNorth(1)], [0 CNorth(2)] );
    LineCSouth = line( [0 CSouth(1)], [0 CSouth(2)] );
    set([LineCNorth, LineCSouth], 'LineStyle', '--', 'LineWidth', 3, 'Color',  horzcat(ColorNMM(16,:), 0.7))
    
    % SubSolar direction
    XSun = 40;
    YSun = 40 * tan( acos(ParaSystem.NoseDirection(1)));
    SunDirection = line( [0 XSun], [0 YSun] );
    set([SunDirection], 'LineStyle', ':', 'LineWidth', 3, 'Color', [0 0 0 0.2])
    
    % Nose
    Nose = scatter(25, 0, 200, 'filled');
    ColorNose = cbrewer2('OrRd', 20);
    set(Nose, 'CData', ColorNose(19, :));

    % Options
        axes = gca;
        set(axes, 'FontName'   , 'Helvetica',...
            'XAxisLocation', 'bottom',...
            'YAxisLocation', 'left',...
            'xdir', 'reverse',...
            'FontSize', 16,...
            'LineWidth'   , 1 ...
            );
        Ylabel = ylabel('Z (R_p)');
        Xlabel = xlabel('X (R_p)');
        ColourBackground = [0, 0, 0, 0.3];
        set(gca,'Color',ColourBackground)
        set(gcf,'color','w');
        
    % Axes
    axis([-5 30 -32 32])
    
    export_fig NMM_Phi0.pdf -opengl -q101
    
    (length(ParaGrid.PhiList)*length(ParaGrid.ThetaList))^2
    