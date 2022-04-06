
function [MagField] = MagField_2(rGrid, ParaGrid, ParaSystem_ini)


    global ParaSystem
    ParaSystem = ParaSystem_ini;

%% Extracting Parameters
    ThetaGrid_ini = ParaGrid.ThetaGrid;
    PhiGrid_ini = ParaGrid.PhiGrid;
    MTilted = ParaSystem.MTilted;
    P_CartNose2KSM = ParaSystem.P_CartNose2KSM;
    
%% Building Entire Surface    
    rWhole = vertcat(rGrid, flipud( rGrid(2:end, :)) );
    PhiListWhole = vertcat(ParaGrid.PhiList, pi+ParaGrid.PhiList(2:end, :) );
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

    YGrid_Norm(1, :) = 0;
    YGrid_Norm(end, :) = 0;
    
        
% Magnetic Pressure: Grid-Computation of the Dipole Field
    MdotEr = MTilted(1) .* XGrid_Norm + MTilted(2) .* YGrid_Norm + MTilted(3) .* ZGrid_Norm;

    BDipole.X = (1./rWhole.^3) .* ( 3*MdotEr .* XGrid_Norm - MGridX );
    BDipole.Y = (1./rWhole.^3) .* ( 3*MdotEr .* YGrid_Norm - MGridY );
    BDipole.Z = (1./rWhole.^3) .* ( 3*MdotEr .* ZGrid_Norm - MGridZ );

    BDipole.XInterp = griddedInterpolant( PhiGrid, ThetaGrid, BDipole.X );
    BDipole.YInterp = griddedInterpolant( PhiGrid, ThetaGrid, BDipole.Y );
    BDipole.ZInterp = griddedInterpolant( PhiGrid, ThetaGrid, BDipole.Z );
    

%% Adding CAN-disk Field

% Over Half MP Surface
BCAN = ParaSystem.BCAN;
Elevation = ( pi/2 - ThetaGrid_ini );
[Y_Nose, Z_Nose, X_Nose] = sph2cart(PhiGrid_ini, Elevation, rGrid .* ParaSystem.r0 ./ ParaSystem.Rp);

% [Y_Nose, Z_Nose, X_Nose] = sph2cart(PhiGrid_ini, pi/2-ThetaGrid_ini, rGrid .* ParaSystem.r0 ./ ParaSystem.Rp);

    X_KSM = P_CartNose2KSM(1,1) .* X_Nose + P_CartNose2KSM(1,2) .* Y_Nose + P_CartNose2KSM(1,3) .* Z_Nose;
    Y_KSM = P_CartNose2KSM(2,1) .* X_Nose + P_CartNose2KSM(2,2) .* Y_Nose + P_CartNose2KSM(2,3) .* Z_Nose;
    Z_KSM = P_CartNose2KSM(3,1) .* X_Nose + P_CartNose2KSM(3,2) .* Y_Nose + P_CartNose2KSM(3,3) .* Z_Nose;

    %{
    z = Z_KSM;
    rho = sqrt( X_KSM.^2 + Y_KSM.^2 );
    x = ParaSystem.CAN_DiskParameters;
    xdata = [z(:).'; rho(:).'];
    
    B_CAN = diskfield_cyl(x, xdata);
    B_CAN_Z = B_CAN(1:length(B_CAN)/2) .';
    B_CAN_rho = B_CAN(length(B_CAN)/2+1:length(B_CAN)) .';
    B_CAN_phi = 0.* B_CAN_Z;
    
    B_CAN_Z( PhiGrid_ini(:) == -pi/2 & ThetaGrid_ini(:) == pi/2 ) = - B_CAN_Z( PhiGrid_ini(:) == pi/2 & ThetaGrid_ini(:) == pi/2 );  
    B_CAN_rho( abs(PhiGrid_ini(:)) == pi/2 ) = 0;
    
    PhiList = linspace(0, 2*pi, length(z(:))) .';
    B_CAN_X = cos(PhiList) .* B_CAN_rho - 0 ;
    B_CAN_Y = sin(PhiList) .* B_CAN_rho + 0 ;
    
    B_CAN_Cartesian.X = reshape(B_CAN_X, size(rGrid));
    B_CAN_Cartesian.Y = reshape(B_CAN_Y, size(rGrid));
    B_CAN_Cartesian.Z = reshape(B_CAN_Z, size(rGrid));
    
    BCAN_X = B_CAN_Cartesian.X;
    BCAN_Y = B_CAN_Cartesian.Y;
    BCAN_Z = B_CAN_Cartesian.Z;
%}
    
BCAN_X = BCAN.XInterp( X_KSM, Y_KSM, Z_KSM );
BCAN_Y = BCAN.YInterp( X_KSM, Y_KSM, Z_KSM );
BCAN_Z = BCAN.ZInterp( X_KSM, Y_KSM, Z_KSM );
% 
% B_DipCAN.X = griddedInterpolant(PhiGrid, ThetaGrid, BDipole.X + BCAN_X);
% B_DipCAN.Y = griddedInterpolant(PhiGrid, ThetaGrid, BDipole.Y + BCAN_Y);
% B_DipCAN.Z = griddedInterpolant(PhiGrid, ThetaGrid, BDipole.Z + BCAN_Z);

    
    %% Induced Shielding Field

%     BTotal = ShieldingField(rGrid, ParaGrid, ParaSystem, BDipole) + MagField(rGrid, ParaGrid, ParaSystem)
%  
%     MagField_k = MagField_km1 + BShielding(MagField_km1)

    C_CAN = ParaSystem.OnOff_CAN;

if ParaSystem.OnOff_ShieldingField == 0
    
    BFieldX = BDipole.XInterp(PhiGrid_ini, ThetaGrid_ini) + C_CAN .* BCAN_X;
    BFieldY = BDipole.YInterp(PhiGrid_ini, ThetaGrid_ini) + C_CAN .* BCAN_Y;
    BFieldZ = BDipole.ZInterp(PhiGrid_ini, ThetaGrid_ini) + C_CAN .* BCAN_Z;

    
else
    
    BShielding = ParaSystem.BShielding;
   
    % CAN Disk over WHOLE MP SUrface    
    Elevation = ( pi/2 - ThetaGrid );
    [Y_Nose_Whole, Z_Nose_Whole, X_Nose_Whole] = sph2cart(PhiGrid, Elevation, rWhole .* ParaSystem.r0 ./ ParaSystem.Rp);

    X_KSM_Whole = P_CartNose2KSM(1,1) .* X_Nose_Whole + P_CartNose2KSM(1,2) .* Y_Nose_Whole + P_CartNose2KSM(1,3) .* Z_Nose_Whole;
    Y_KSM_Whole = P_CartNose2KSM(2,1) .* X_Nose_Whole + P_CartNose2KSM(2,2) .* Y_Nose_Whole + P_CartNose2KSM(2,3) .* Z_Nose_Whole;
    Z_KSM_Whole = P_CartNose2KSM(3,1) .* X_Nose_Whole + P_CartNose2KSM(3,2) .* Y_Nose_Whole + P_CartNose2KSM(3,3) .* Z_Nose_Whole;

    B_CAN_Whole.X = BCAN.XInterp( X_KSM_Whole, Y_KSM_Whole, Z_KSM_Whole );
    B_CAN_Whole.Y = BCAN.XInterp( X_KSM_Whole, Y_KSM_Whole, Z_KSM_Whole );
    B_CAN_Whole.Z = BCAN.XInterp( X_KSM_Whole, Y_KSM_Whole, Z_KSM_Whole );

    Method = 'spline';
    B_DipCAN.XInterp = griddedInterpolant(PhiGrid, ThetaGrid, BDipole.X + C_CAN .* B_CAN_Whole.X, Method);
    B_DipCAN.YInterp = griddedInterpolant(PhiGrid, ThetaGrid, BDipole.Y + C_CAN .* B_CAN_Whole.Y, Method);
    B_DipCAN.ZInterp = griddedInterpolant(PhiGrid, ThetaGrid, BDipole.Z + C_CAN .* B_CAN_Whole.Z, Method);

        if not( isequal([BShielding.X; BShielding.Y; BShielding.Z], [0; 0; 0]) )

    %         BShielding.X = BShielding.XInterp_r(rGrid);
    %         BShielding.Y = BShielding.YInterp_r(rGrid);
    %         BShielding.Z = BShielding.ZInterp_r(rGrid);

    %         BShielding.X = BShielding.XInterp(PhiGrid, ThetaGrid);
    %         BShielding.Y = BShielding.YInterp(PhiGrid, ThetaGrid);
    %         BShielding.Z = BShielding.ZInterp(PhiGrid, ThetaGrid);
            

    %         BTotal.X = BDipole.X + BShielding.XInterp(PhiGrid, ThetaGrid);
            BTotal.X = BDipole.X + C_CAN .* B_CAN_Whole.X + BShielding.XInterp(PhiGrid, ThetaGrid);
            BTotal.XInterp = griddedInterpolant(PhiGrid, ThetaGrid, BTotal.X);

    %         BTotal.Y = BDipole.Y + BShielding.YInterp(PhiGrid, ThetaGrid);
            BTotal.Y = BDipole.Y + C_CAN .* B_CAN_Whole.Y + BShielding.YInterp(PhiGrid, ThetaGrid);
            BTotal.YInterp = griddedInterpolant(PhiGrid, ThetaGrid, BTotal.Y);

    %         BTotal.Z = BDipole.Z + BShielding.ZInterp(PhiGrid, ThetaGrid);
            BTotal.Z = BDipole.Z + C_CAN .* B_CAN_Whole.Z + BShielding.ZInterp(PhiGrid, ThetaGrid);
            BTotal.ZInterp = griddedInterpolant(PhiGrid, ThetaGrid, BTotal.Z);

%             BShielding = ShieldingField(rGrid, ParaGrid, ParaSystem, BTotal);
%             ParaSystem.BShielding = BShielding;

ParaSystem.BShielding = ShieldingField(rGrid, ParaGrid, ParaSystem, BTotal);


        else

    %        BShielding = ShieldingField(rGrid, ParaGrid, ParaSystem, BDipole);
           BShielding = ShieldingField(rGrid, ParaGrid, ParaSystem, B_DipCAN);
           ParaSystem.BShielding = BShielding;

        end

    BFieldX = BDipole.XInterp(PhiGrid_ini, ThetaGrid_ini) + C_CAN .* BCAN_X + BShielding.XInterp(PhiGrid_ini, ThetaGrid_ini);
    BFieldY = BDipole.YInterp(PhiGrid_ini, ThetaGrid_ini) + C_CAN .* BCAN_Y + BShielding.YInterp(PhiGrid_ini, ThetaGrid_ini);
    BFieldZ = BDipole.ZInterp(PhiGrid_ini, ThetaGrid_ini) + C_CAN .* BCAN_Z + BShielding.ZInterp(PhiGrid_ini, ThetaGrid_ini);

end

%     BFieldX = BDipole.XInterp(PhiGrid_ini, ThetaGrid_ini) + BShielding.XInterp(PhiGrid_ini, ThetaGrid_ini);
%     BFieldY = BDipole.YInterp(PhiGrid_ini, ThetaGrid_ini) + BShielding.YInterp(PhiGrid_ini, ThetaGrid_ini);
%     BFieldZ = BDipole.ZInterp(PhiGrid_ini, ThetaGrid_ini) + BShielding.ZInterp(PhiGrid_ini, ThetaGrid_ini);
    %}
    
%     


% %     
% 
   %{
    rWhole(2, :) = rWhole(end-1, :);
    rWhole(1, :) = (1/2) .* ( rWhole(2, :) + rWhole(end-1, :) );
    rWhole(end, :) = rWhole(1, :);
    [Y, Z, X] = sph2cart(PhiGrid, pi/2 - ThetaGrid, rWhole .* ParaSystem.r0 ./ ParaSystem.Rp);

    BWhole_X = BTotal.XInterp(PhiGrid, ThetaGrid);
    BWhole_Y = BTotal.YInterp(PhiGrid, ThetaGrid);
    BWhole_Z = BTotal.ZInterp(PhiGrid, ThetaGrid);
    
    BDip_X = BDipole.XInterp(PhiGrid, ThetaGrid);
    BDip_Y = BDipole.YInterp(PhiGrid, ThetaGrid);
    BDip_Z = BDipole.ZInterp(PhiGrid, ThetaGrid);

    XKSM = P_CartNose2KSM(1,1) .* X + P_CartNose2KSM(1,2) .* Y + P_CartNose2KSM(1,3) .* Z;
    YKSM = P_CartNose2KSM(2,1) .* X + P_CartNose2KSM(2,2) .* Y + P_CartNose2KSM(2,3) .* Z;
    ZKSM = P_CartNose2KSM(3,1) .* X + P_CartNose2KSM(3,2) .* Y + P_CartNose2KSM(3,3) .* Z;

    figure;
    C = 1;
    hold on
        Quiver1 = quiver3(XKSM(1:C:end,1:C:end) , YKSM(1:C:end,1:C:end), ZKSM(1:C:end,1:C:end), BWhole_X(1:C:end,1:C:end), BWhole_Y(1:C:end,1:C:end), BWhole_Z(1:C:end,1:C:end));
        Quiver2 = quiver3(XKSM, YKSM, ZKSM, BDip_X, BDip_Y, BDip_Z);
        set([Quiver1, Quiver2], 'AutoScaleFactor', 2, 'LineWidth', 3, 'MaxHeadSize', 1)
        axis equal
        axis([-5 30 -25 25 -23 23]);
        PlotSurf = surf(XKSM, YKSM, ZKSM, XKSM);
        set(PlotSurf, 'FaceAlpha', 0.3)
    %}
    


%% Total Field
% C = ParaSystem.OnOff_CAN;
%     BFieldX = BDipole.XInterp(PhiGrid_ini, ThetaGrid_ini) + C .* BCAN_X;
%     BFieldY = BDipole.YInterp(PhiGrid_ini, ThetaGrid_ini) + C .* BCAN_Y;
%     BFieldZ = BDipole.ZInterp(PhiGrid_ini, ThetaGrid_ini) + C .* BCAN_Z;

% 
%     BFieldX = BDipole.XInterp(PhiGrid_ini, ThetaGrid_ini) + C .* BCAN_X + BShielding.XInterp(PhiGrid_ini, ThetaGrid_ini);
%     BFieldY = BDipole.YInterp(PhiGrid_ini, ThetaGrid_ini) + C .* BCAN_Y + BShielding.YInterp(PhiGrid_ini, ThetaGrid_ini);
%     BFieldZ = BDipole.ZInterp(PhiGrid_ini, ThetaGrid_ini) + C .* BCAN_Z + BShielding.ZInterp(PhiGrid_ini, ThetaGrid_ini);
%}
    
% 
% XPlot = ThetaGrid_ini;
% YPlot = PhiGrid_ini;
% ZPlot = BDipole.ZInterp(PhiGrid_ini, ThetaGrid_ini) ;
% scatter3(XPlot(:), YPlot(:), ZPlot(:))

% % 
% figure;
%     pcolor((BFieldZ));
% colorbar;

%     BFieldX = BDipole.XInterp(PhiGrid_ini, ThetaGrid_ini);
%     BFieldY = BDipole.YInterp(PhiGrid_ini, ThetaGrid_ini);
%     BFieldZ = BDipole.ZInterp(PhiGrid_ini, ThetaGrid_ini);
    

%% Packaging Output

MagField.X = BFieldX;
MagField.Y = BFieldY;
MagField.Z = BFieldZ;

MagField.XInterp = griddedInterpolant(PhiGrid_ini, ThetaGrid_ini, BFieldX);
MagField.YInterp = griddedInterpolant(PhiGrid_ini, ThetaGrid_ini, BFieldY);
MagField.ZInterp = griddedInterpolant(PhiGrid_ini, ThetaGrid_ini, BFieldZ);


end