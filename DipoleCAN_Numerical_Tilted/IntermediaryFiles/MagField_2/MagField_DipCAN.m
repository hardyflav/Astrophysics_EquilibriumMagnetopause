

function [B_DipCAN] = MagField_DipCAN(rGrid, ParaGrid, ParaSystem_ini)

    global ParaSystem
    ParaSystem = ParaSystem_ini;

%% Extracting Parameters
    ThetaGrid_ini = ParaGrid.ThetaGrid;
    PhiGrid_ini = ParaGrid.PhiGrid;
    MTilted = ParaSystem.MTilted;
    P_CartNose2KSM = ParaSystem.P_CartNose2KSM;
    
    
%% Building Entire Surface    
    rWhole = vertcat(rGrid, flipud( rGrid(1:end-1, :)) );
%     PhiListWhole = vertcat(ParaGrid.PhiList, pi+ParaGrid.PhiList(2:end-1, :) );
    PhiListWhole = linspace(-pi/2, 3*pi/2, size(rWhole, 1));
    [PhiGrid, ThetaGrid] = ndgrid( PhiListWhole, ParaGrid.ThetaList );
    
    
%% Dipole Field

% Magnetic Pressure: Rotated Magnetic Moment, KSM
    MGridX = 0*rWhole + MTilted(1);
    MGridY = 0*rWhole + MTilted(2);
    MGridZ = 0*rWhole + MTilted(3);

% Magnetic Pressure: Position of Nose, KSM
    XGrid_Norm_Nose = 1 .* cos(ThetaGrid);
    YGrid_Norm_Nose = 1 .* sin(ThetaGrid) .* cos(PhiGrid);
    ZGrid_Norm_Nose = 1 .* sin(ThetaGrid) .* sin(PhiGrid);
    
    XGrid_Norm = P_CartNose2KSM(1,1) .* XGrid_Norm_Nose + P_CartNose2KSM(1,2) .* YGrid_Norm_Nose + P_CartNose2KSM(1,3) .* ZGrid_Norm_Nose;
    YGrid_Norm = P_CartNose2KSM(2,1) .* XGrid_Norm_Nose + P_CartNose2KSM(2,2) .* YGrid_Norm_Nose + P_CartNose2KSM(2,3) .* ZGrid_Norm_Nose;
    ZGrid_Norm = P_CartNose2KSM(3,1) .* XGrid_Norm_Nose + P_CartNose2KSM(3,2) .* YGrid_Norm_Nose + P_CartNose2KSM(3,3) .* ZGrid_Norm_Nose;

%     YGrid_Norm(1, :) = 0;
%     YGrid_Norm(end, :) = 0;
    
        
% Magnetic Pressure: Grid-Computation of the Dipole Field
    MdotEr = MTilted(1) .* XGrid_Norm + MTilted(2) .* YGrid_Norm + MTilted(3) .* ZGrid_Norm;

    BDipole.X = (1./rWhole.^3) .* ( 3*MdotEr .* XGrid_Norm - MGridX );
    BDipole.Y = (1./rWhole.^3) .* ( 3*MdotEr .* YGrid_Norm - MGridY );
    BDipole.Z = (1./rWhole.^3) .* ( 3*MdotEr .* ZGrid_Norm - MGridZ );

    %{
    BDipole.XInterp = griddedInterpolant( PhiGrid, ThetaGrid, BDipole.X );
    BDipole.YInterp = griddedInterpolant( PhiGrid, ThetaGrid, BDipole.Y );
    BDipole.ZInterp = griddedInterpolant( PhiGrid, ThetaGrid, BDipole.Z );
    %}

    
%% Adding CAN-disk Field

% Over Entire MP Surface
    BCAN = ParaSystem.BCAN;
    Elevation = ( pi/2 - ThetaGrid );
    [Y_Nose, Z_Nose, X_Nose] = sph2cart(PhiGrid, Elevation, rWhole .* ParaSystem.r0 ./ ParaSystem.Rp);
    X_KSM = P_CartNose2KSM(1,1) .* X_Nose + P_CartNose2KSM(1,2) .* Y_Nose + P_CartNose2KSM(1,3) .* Z_Nose;
    Y_KSM = P_CartNose2KSM(2,1) .* X_Nose + P_CartNose2KSM(2,2) .* Y_Nose + P_CartNose2KSM(2,3) .* Z_Nose;
    Z_KSM = P_CartNose2KSM(3,1) .* X_Nose + P_CartNose2KSM(3,2) .* Y_Nose + P_CartNose2KSM(3,3) .* Z_Nose;

    BCAN_X = BCAN.XInterp( X_KSM, Y_KSM, Z_KSM );
    BCAN_Y = BCAN.YInterp( X_KSM, Y_KSM, Z_KSM );
    BCAN_Z = BCAN.ZInterp( X_KSM, Y_KSM, Z_KSM );

    Method = 'linear';
    B_DipCAN.XInterp = griddedInterpolant(PhiGrid, ThetaGrid, BDipole.X + BCAN_X, Method);
    B_DipCAN.YInterp = griddedInterpolant(PhiGrid, ThetaGrid, BDipole.Y + BCAN_Y, Method);
    B_DipCAN.ZInterp = griddedInterpolant(PhiGrid, ThetaGrid, BDipole.Z + BCAN_Z, Method);


end