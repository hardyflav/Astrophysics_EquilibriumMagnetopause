
function [MagField] = MagField(rGrid, ParaGrid, ParaSystem_ini)


    global ParaSystem
    ParaSystem = ParaSystem_ini;

%% Extracting Parameters
    ThetaGrid_ini = ParaGrid.ThetaGrid;
    PhiGrid_ini = ParaGrid.PhiGrid;
    MTilted = ParaSystem.MTilted;
    P_CartNose2KSM = ParaSystem.P_CartNose2KSM;
    
%% Building Entire Surface    
%     rWhole = vertcat(rGrid, flipud( rGrid(1:end-1, :)) );
%     PhiListWhole = vertcat(ParaGrid.PhiList, pi+ParaGrid.PhiList(1:end-1, :) );
    rWhole = vertcat(rGrid, flipud( rGrid(1:end-1, :)) );
%     PhiListWhole = vertcat(ParaGrid.PhiList, pi+ParaGrid.PhiList(2:end-1, :) );
    PhiListWhole = linspace(-pi/2, 3*pi/2, size(rWhole, 1));

    [PhiGrid, ThetaGrid] = ndgrid( PhiListWhole, ParaGrid.ThetaList );
    
    Elevation = ( pi/2 - ThetaGrid );
    [Y_Nose, Z_Nose, X_Nose] = sph2cart(PhiGrid, Elevation, rWhole .* ParaSystem.r0 ./ ParaSystem.Rp);
    X_KSM = P_CartNose2KSM(1,1) .* X_Nose + P_CartNose2KSM(1,2) .* Y_Nose + P_CartNose2KSM(1,3) .* Z_Nose;
    Y_KSM = P_CartNose2KSM(2,1) .* X_Nose + P_CartNose2KSM(2,2) .* Y_Nose + P_CartNose2KSM(2,3) .* Z_Nose;
    Z_KSM = P_CartNose2KSM(3,1) .* X_Nose + P_CartNose2KSM(3,2) .* Y_Nose + P_CartNose2KSM(3,3) .* Z_Nose;
    
    
%% Dip + CAN Field over WHOLE surface
    B_DipCAN = MagField_DipCAN(rGrid, ParaGrid, ParaSystem) ;   
    B_DipCAN.XValWhole = B_DipCAN.XInterp(PhiGrid, ThetaGrid);
    B_DipCAN.YValWhole = B_DipCAN.YInterp(PhiGrid, ThetaGrid);
    B_DipCAN.ZValWhole = B_DipCAN.ZInterp(PhiGrid, ThetaGrid);

    B_DipCAN.XHalf = B_DipCAN.XInterp(ParaGrid.PhiGrid, ParaGrid.ThetaGrid);
    B_DipCAN.YHalf = B_DipCAN.YInterp(ParaGrid.PhiGrid, ParaGrid.ThetaGrid);
    B_DipCAN.ZHalf = B_DipCAN.ZInterp(ParaGrid.PhiGrid, ParaGrid.ThetaGrid);
    
       %{
       % Visualisation
       figure;
       quiver3(X_KSM, Y_KSM, Z_KSM, B_DipCAN.XValWhole, B_DipCAN.YValWhole, B_DipCAN.ZValWhole);
       axis equal
        %}

%% Shielding Field over WHOLE surface

%{
% TEST: INTERPOLATION
    if ParaSystem.OnOff_ShieldingField == 1
        BShielding = ParaSystem.BShielding;
        BShielding.XValWhole = BShielding.XInterp(X_KSM, Y_KSM, Z_KSM);
        BShielding.YValWhole = BShielding.YInterp(X_KSM, Y_KSM, Z_KSM);
        BShielding.ZValWhole = BShielding.ZInterp(X_KSM, Y_KSM, Z_KSM);
        
        BTot.XValWhole = B_DipCAN.XValWhole + BShielding.XValWhole;
        BTot.YValWhole = B_DipCAN.YValWhole + BShielding.YValWhole;
        BTot.ZValWhole = B_DipCAN.ZValWhole + BShielding.ZValWhole;
    else
        BTot.XValWhole = B_DipCAN.XValWhole;
        BTot.YValWhole = B_DipCAN.YValWhole;
        BTot.ZValWhole = B_DipCAN.ZValWhole;
    end
%}    

if  ParaSystem.OnOff_ShieldingField == 0

    BTot.X = B_DipCAN.XHalf;
    BTot.Y = B_DipCAN.YHalf;
    BTot.Z = B_DipCAN.ZHalf;

else

% Building Half Surface
%{
    Elevation = ( pi/2 - ParaGrid.ThetaGrid );
    [Y_Nose, Z_Nose, X_Nose] = sph2cart(ParaGrid.PhiGrid, Elevation, rGrid .* ParaSystem.r0 ./ ParaSystem.Rp);
    X_KSM = P_CartNose2KSM(1,1) .* X_Nose + P_CartNose2KSM(1,2) .* Y_Nose + P_CartNose2KSM(1,3) .* Z_Nose;
    Y_KSM = P_CartNose2KSM(2,1) .* X_Nose + P_CartNose2KSM(2,2) .* Y_Nose + P_CartNose2KSM(2,3) .* Z_Nose;
    Z_KSM = P_CartNose2KSM(3,1) .* X_Nose + P_CartNose2KSM(3,2) .* Y_Nose + P_CartNose2KSM(3,3) .* Z_Nose;
%}
    CPhi = 1;
    CTheta = 1;
    MListX_temp = X_KSM(1:CPhi:end, 1:CTheta:end);
    MList.X = MListX_temp(:);
    MListY_temp = Y_KSM(1:CPhi:end, 1:CTheta:end);
    MList.Y = MListY_temp(:);
    MListZ_temp = Z_KSM(1:CPhi:end, 1:CTheta:end);
    MList.Z = MListZ_temp(:);

    BM_X_List = 0.*MList.X;
    BM_Y_List = 0.*BM_X_List;
    BM_Z_List = 0.*BM_X_List;
    
    if isequal([ParaSystem.BShielding.XVal, ParaSystem.BShielding.YVal, ParaSystem.BShielding.ZVal], [0, 0, 0])
        BField_km1 = B_DipCAN ;
    else
        BShield_km1 = ParaSystem.BShielding ;
        BField_km1.XValWhole = B_DipCAN.XValWhole + BShield_km1.XValWhole;
        BField_km1.YValWhole = B_DipCAN.YValWhole + BShield_km1.YValWhole;
        BField_km1.ZValWhole = B_DipCAN.ZValWhole + BShield_km1.ZValWhole;
        Method = 'linear';
        BField_km1.XInterp = griddedInterpolant(PhiGrid, ThetaGrid, BField_km1.XValWhole, Method);
        BField_km1.YInterp = griddedInterpolant(PhiGrid, ThetaGrid, BField_km1.YValWhole, Method);
        BField_km1.ZInterp = griddedInterpolant(PhiGrid, ThetaGrid, BField_km1.ZValWhole, Method);
    end
    
    
      
    for kMx = 1:length(MList.X)
        100*(kMx/length(MList.X));
        M.X = MList.X(kMx);
        M.Y = MList.Y(kMx);
        M.Z = MList.Z(kMx);
        
%         BShielding = ParaSystem.BShielding;
%         BM_X_List(kMx) = BShielding.XInterp( M.X, M.Y, M.Z );
%         BM_Y_List(kMx) = BShielding.YInterp( M.X, M.Y, M.Z );
%         BM_Z_List(kMx) = BShielding.ZInterp( M.X, M.Y, M.Z );
        
        PMx = M.X - X_KSM ;
        PMy = M.Y - Y_KSM ;
        PMz = M.Z - Z_KSM ;
        PM = sqrt(PMx.^2+PMy.^2+PMz.^2);
        
        Cond = ( PM > -1);
        
 
        
        dBpX = @(X) dBpM(rGrid, [X(:,1), X(:,2)], [M.X, M.Y, M.Z], ParaGrid, ParaSystem_ini, BField_km1) ;
        
        PhiP = PhiGrid(Cond);
        ThetaP = ThetaGrid(Cond);

        dBp = dBpX( [PhiP(:), ThetaP(:)] );
        
        
        %{
      figure;   
      hold on
      Plot = scatter3(PhiP, ThetaP, dBp(:,3), 100, PM(Cond), 'filled');
      axis square
      colorbar
        %}
        
%         BM_X_List(kMx) = sum( dBp(:,1) );
%         BM_Y_List(kMx) = sum( dBp(:,2) );
%         BM_Z_List(kMx) = sum( dBp(:,3) );
        
        %{
        % TEST: Monte-Carlo Integrators
        dBpX = @(X) dBpM(rGrid, [X(:,1), X(:,2)], [M.X, M.Y, M.Z], ParaGrid, ParaSystem_ini, BField_km1, 'X') ;
        dBpY = @(X) dBpM(rGrid, [X(:,1), X(:,2)], [M.X, M.Y, M.Z], ParaGrid, ParaSystem_ini, BField_km1, 'Y') ;
        dBpZ = @(X) dBpM(rGrid, [X(:,1), X(:,2)], [M.X, M.Y, M.Z], ParaGrid, ParaSystem_ini, BField_km1, 'Z') ;
        BM_X_List(kMx) = integralN_mc(dBpX, [-pi/2, 3*pi/2; 0, ParaGrid.ThetaList(end)], 'timelimit', 0.1) ;
        BM_Y_List(kMx) = integralN_mc(dBpY, [-pi/2, 3*pi/2; 0, ParaGrid.ThetaList(end)], 'timelimit', 0.1) ;
        BM_Z_List(kMx) = integralN_mc(dBpZ, [-pi/2, 3*pi/2; 0, ParaGrid.ThetaList(end)], 'timelimit', 0.1) ;
        %}
        
% ------ TEST:  Trapezoidal Rule
        PhiVal_Interp = PhiGrid(1:CPhi:end, 1:CTheta:end);
        ThetaVal_Interp = ThetaGrid(1:CPhi:end, 1:CTheta:end);
        
        dBpX_Grid = reshape( dBp(:, 1), size(PhiVal_Interp) );
        dBpY_Grid = reshape( dBp(:, 2), size(PhiVal_Interp) );
        dBpZ_Grid = reshape( dBp(:, 3), size(PhiVal_Interp) );
        
        %{
        Qx_Phi = dBpX_Grid(1, :) ./ 2 + sum(dBpX_Grid(2:end-1, :)) + dBpX_Grid(end, :) ./ 2;
        Qx_Theta = Qx_Phi(1)/2 + sum(Qx_Phi(2:end-1)) + Qx_Phi(end)/2;
        BM_X_List(kMx) = Qx_Theta;
        
        Qy_Phi = dBpY_Grid(1, :) ./ 2 + sum(dBpY_Grid(2:end-1, :)) + dBpY_Grid(end, :) ./ 2;
        Qy_Theta = Qy_Phi(1)/2 + sum(Qy_Phi(2:end-1)) + Qy_Phi(end)/2;
        BM_Y_List(kMx) = Qy_Theta;
        
        Qz_Phi = dBpZ_Grid(1, :) ./ 2 + sum(dBpZ_Grid(2:end-1, :)) + dBpZ_Grid(end, :) ./ 2;
        Qz_Theta = Qz_Phi(1)/2 + sum(Qz_Phi(2:end-1)) + Qz_Phi(end)/2;
        BM_Z_List(kMx) = Qz_Theta;
        %}
        
        
        
        
        
        
        BM_X_List(kMx) = (1/4) .* ( dBpX_Grid(1, 1) + dBpX_Grid(end, 1) + dBpX_Grid(1, end) + dBpX_Grid(end, end) + ...
                           2 .* sum(dBpX_Grid(2:end-1, 1)) + 2 .* sum(dBpX_Grid(2:end-1, end)) + 2 .* sum(dBpX_Grid(1, 2:end-1)) + 2 .* sum(dBpX_Grid(end, 2:end-1)) + ...
                           4 .* sum( sum(dBpX_Grid(2:end-1, 2:end-1), 1)) );
        
        BM_Y_List(kMx) = (1/4) .* ( dBpY_Grid(1, 1) + dBpY_Grid(end, 1) + dBpY_Grid(1, end) + dBpY_Grid(end, end) + ...
                           2 .* sum(dBpY_Grid(2:end-1, 1)) + 2 .* sum(dBpY_Grid(2:end-1, end)) + 2 .* sum(dBpY_Grid(1, 2:end-1)) + 2 .* sum(dBpY_Grid(end, 2:end-1)) + ...
                           4 .* sum( sum(dBpY_Grid(2:end-1, 2:end-1), 1)) );
               
        BM_Z_List(kMx) = (1/4) .* ( dBpZ_Grid(1, 1) + dBpZ_Grid(end, 1) + dBpZ_Grid(1, end) + dBpZ_Grid(end, end) + ...
                           2 .* sum(dBpZ_Grid(2:end-1, 1)) + 2 .* sum(dBpZ_Grid(2:end-1, end)) + 2 .* sum(dBpZ_Grid(1, 2:end-1)) + 2 .* sum(dBpZ_Grid(end, 2:end-1)) + ...
                           4 .* sum( sum(dBpZ_Grid(2:end-1, 2:end-1), 1) ) );
               
                       
        % TEST: INTERPOLATION
        %{
        dBpZ_Grid2 = dBpZ_Grid;
        Cond = PM < 20;
        dBpZ_Grid = dBpZ_Grid2 .* Cond;
        
        
        [row, col] = find( PM == 0 );
        
        Method = 'linear';
        dBpX_Interp = griddedInterpolant( PhiGrid, ThetaGrid, dBpX_Grid, Method );
        dBpY_Interp = griddedInterpolant( PhiGrid, ThetaGrid, dBpY_Grid, Method );
        dBpZ_Interp = griddedInterpolant( PhiGrid, ThetaGrid, dBpZ_Grid, Method );
        PM_Interp = griddedInterpolant(PhiGrid, ThetaGrid, PM, Method);
       
        hold on
       figure;
       scatter3(PhiGrid(:), ThetaGrid(:), dBpZ_Grid(:), 40, PM(:), 'filled')
        colorbar
        
        
        
        
       PhiList_Interp = linspace(PhiListWhole(1), PhiListWhole(end), length(ParaGrid.PhiList) .* 1);
       ThetaList_Interp = linspace(ParaGrid.ThetaList(1), ParaGrid.ThetaList(end), length(ParaGrid.ThetaList) .* 1);
       [PhiGrid_Interp, ThetaGrid_Interp] = ndgrid(PhiList_Interp, ThetaList_Interp);
        
       PM_Val = PM_Interp(PhiGrid_Interp, ThetaGrid_Interp);
       dBpZ_Val = dBpZ_Interp(PhiGrid_Interp, ThetaGrid_Interp);
        [row, col] = find( PM_Val == 0 );
        
        
       scatter3(PhiGrid_Interp(:), ThetaGrid_Interp(:), dBpZ_Val(:), 'filled')
        
        
        
        
(1/4) .* ( dBpZ_Val(1, 1) + dBpZ_Val(end, 1) + dBpZ_Val(1, end) + dBpZ_Val(end, end) + ...
                           2 .* sum(dBpZ_Val(2:end-1, 1)) + 2 .* sum(dBpZ_Val(2:end-1, end)) + 2 .* sum(dBpZ_Val(1, 2:end-1)) + 2 .* sum(dBpZ_Val(end, 2:end-1)) + ...
                           4 .* sum( sum(dBpZ_Val(2:end-1, 2:end-1), 1) ) )
        
        
        figure;
        pcolor(dBpZ_Val); colorbar
        shading flat        
        
        
        Integ = @(Phi, Theta) dBpx_Interp(Phi, Theta);
        integral2(Integ, -pi/2, 3*pi/2, )
        
        
        
       Cond = abs(dBpZ_Grid) > max( abs(dBpZ_Grid(:)) ) / 5;
       Cond = PM > 0.1;
       dBpZ_Cont = dBpZ_Grid .* Cond;
       dBpZ_Cont>0
       dBpZ_Cont_Interp = griddedInterpolant(PhiGrid(:, 2:end), ThetaGrid(:, 2:end), dBpZ_Grid(:, 2:end) );
       
       
        (1/4) .* ( dBpZ_Cont_Val(1, 1) + dBpZ_Cont_Val(end, 1) + dBpZ_Cont_Val(1, end) + dBpZ_Cont_Val(end, end) + ...
                           2 .* sum(dBpZ_Cont_Val(2:end-1, 1)) + 2 .* sum(dBpZ_Cont_Val(2:end-1, end)) + 2 .* sum(dBpZ_Cont_Val(1, 2:end-1)) + 2 .* sum(dBpZ_Cont_Val(end, 2:end-1)) + ...
                           4 .* sum( sum(dBpZ_Cont_Val(2:end-1, 2:end-1), 1)) )
                       
        figure;
        pcolor( dBpZ_Grid )
        colorbar               
                       
        trapz( trapz(dBpZ_Cont_Val,1), 2)        
       
       
       figure;
        pcolor(BM_X_Val)               
        colorbar           
        %}
                       
                       
% ------ TEST:  Trapezoidal Rule
        
    end
       %{
       % Visualisation
       figure;
        hold on
       quiver3(X_KSM, Y_KSM, Z_KSM, reshape(BM_X_List, size(PhiVal_Interp) ), reshape(BM_Y_List, size(PhiVal_Interp) ), reshape(BM_Z_List, size(PhiVal_Interp) ));
       quiver3(X_KSM, Y_KSM, Z_KSM, reshape(BM_X_List, size(PhiVal_Interp) ) + B_DipCAN.XValWhole, reshape(BM_Y_List, size(PhiVal_Interp) ) + B_DipCAN.YValWhole, reshape(BM_Z_List, size(PhiVal_Interp) ) + B_DipCAN.ZValWhole);
       axis equal
        %}
    
    
    
    PhiVal_Interp = PhiGrid(1:CPhi:end, 1:CTheta:end);
    ThetaVal_Interp = ThetaGrid(1:CPhi:end, 1:CTheta:end);
    Method = 'linear';
    BShield.XInterp = griddedInterpolant( PhiVal_Interp, ThetaVal_Interp, reshape(BM_X_List, size(PhiVal_Interp) ), Method );
    BShield.YInterp = griddedInterpolant( PhiVal_Interp, ThetaVal_Interp, reshape(BM_Y_List, size(PhiVal_Interp) ), Method );
    BShield.ZInterp = griddedInterpolant( PhiVal_Interp, ThetaVal_Interp, reshape(BM_Z_List, size(PhiVal_Interp) ), Method );

    BShield.XVal = BShield.XInterp( ParaGrid.PhiGrid, ParaGrid.ThetaGrid );
    BShield.YVal = BShield.YInterp( ParaGrid.PhiGrid, ParaGrid.ThetaGrid );
    BShield.ZVal = BShield.ZInterp( ParaGrid.PhiGrid, ParaGrid.ThetaGrid );
    
    BShield.XValWhole = reshape( BM_X_List, size(PhiVal_Interp) );
    BShield.YValWhole = reshape( BM_Y_List, size(PhiVal_Interp) );
    BShield.ZValWhole = reshape( BM_Z_List, size(PhiVal_Interp) );
    
    
    
%{
% Visualisation
    figure;
    pcolor(reshape(BM_Z_List, size(PhiVal_Interp) ))
    pcolor(BShield.ZVal)
    colorbar
    caxis([-1.3 0.7])
%}

    BTot.X = B_DipCAN.XHalf + BShield.XVal;
    BTot.Y = B_DipCAN.YHalf + BShield.YVal;
    BTot.Z = B_DipCAN.ZHalf + BShield.ZVal;
    
    %{
    % Visualisation
        Elevation = ( pi/2 - ParaGrid.ThetaGrid );
        [Y_Nose, Z_Nose, X_Nose] = sph2cart(ParaGrid.PhiGrid, Elevation, rGrid .* ParaSystem.r0 ./ ParaSystem.Rp);
        X_KSM = P_CartNose2KSM(1,1) .* X_Nose + P_CartNose2KSM(1,2) .* Y_Nose + P_CartNose2KSM(1,3) .* Z_Nose;
        Y_KSM = P_CartNose2KSM(2,1) .* X_Nose + P_CartNose2KSM(2,2) .* Y_Nose + P_CartNose2KSM(2,3) .* Z_Nose;
        Z_KSM = P_CartNose2KSM(3,1) .* X_Nose + P_CartNose2KSM(3,2) .* Y_Nose + P_CartNose2KSM(3,3) .* Z_Nose;
    
        figure;
        hold on
        quiver3(X_KSM, Y_KSM, Z_KSM, BTot.X, BTot.Y, BTot.Z)
        quiver3(X_KSM, Y_KSM, Z_KSM, B_DipCAN.XHalf, B_DipCAN.YHalf, B_DipCAN.ZHalf)
    axis equal
    %}

end


%% TOTAL Field over WHOLE surface
    
%{
% TEST: INTERPOLATION
    Method = 'linear';
    BTot.XInterp = griddedInterpolant(PhiGrid, ThetaGrid, BTot.XValWhole, Method);
    BTot.YInterp = griddedInterpolant(PhiGrid, ThetaGrid, BTot.YValWhole, Method);
    BTot.ZInterp = griddedInterpolant(PhiGrid, ThetaGrid, BTot.ZValWhole, Method);

    
    figure;
    hold on
    quiver3(X_KSM, Y_KSM, Z_KSM, BShielding.XValWhole, BShielding.YValWhole, BShielding.ZValWhole)
    axis equal
    
    figure;
    pcolor( Bx )
    caxis([0 0.5])
    colorbar
%}

    Method = 'linear';
    BTot.XInterp = griddedInterpolant( ParaGrid.PhiGrid, ParaGrid.ThetaGrid, BTot.X, Method );
    BTot.YInterp = griddedInterpolant( ParaGrid.PhiGrid, ParaGrid.ThetaGrid, BTot.Y, Method );
    BTot.ZInterp = griddedInterpolant( ParaGrid.PhiGrid, ParaGrid.ThetaGrid, BTot.Z, Method );


%{
% Visualisation

Elevation = ( pi/2 - ParaGrid.ThetaGrid );
[Y_Nose, Z_Nose, X_Nose] = sph2cart(ParaGrid.PhiGrid, Elevation, rGrid .* ParaSystem.r0 ./ ParaSystem.Rp);
X_KSM = P_CartNose2KSM(1,1) .* X_Nose + P_CartNose2KSM(1,2) .* Y_Nose + P_CartNose2KSM(1,3) .* Z_Nose;
Y_KSM = P_CartNose2KSM(2,1) .* X_Nose + P_CartNose2KSM(2,2) .* Y_Nose + P_CartNose2KSM(2,3) .* Z_Nose;
Z_KSM = P_CartNose2KSM(3,1) .* X_Nose + P_CartNose2KSM(3,2) .* Y_Nose + P_CartNose2KSM(3,3) .* Z_Nose;

figure;
hold on
quiver3( X_KSM, Y_KSM, Z_KSM, BTot.X, BTot.Y, BTot.Z )
quiver3( X_KSM, Y_KSM, Z_KSM, B_DipCAN.XHalf, B_DipCAN.YHalf, B_DipCAN.ZHalf )
axis equal
%}

%% Packaging Output
BTot.X = BTot.XInterp( ParaGrid.PhiGrid, ParaGrid.ThetaGrid );
BTot.Y = BTot.YInterp( ParaGrid.PhiGrid, ParaGrid.ThetaGrid );
BTot.Z = BTot.ZInterp( ParaGrid.PhiGrid, ParaGrid.ThetaGrid );

MagField = BTot;

 if  ParaSystem.OnOff_ShieldingField == 1

    ParaSystem.Bshielding_kp1 = BShield;
    ParaSystem.BShielding = ParaSystem.Bshielding_kp1;

 end
    
    
end    
    
    
    %% SIDE
    %{
    %{
 %{    
    
    quiver3(X_KSM, Y_KSM, Z_KSM, )
    
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

%             BShielding.X = BShielding.XInterp_r(rGrid);
%             BShielding.Y = BShielding.YInterp_r(rGrid);
%             BShielding.Z = BShielding.ZInterp_r(rGrid);

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



    
    %}
