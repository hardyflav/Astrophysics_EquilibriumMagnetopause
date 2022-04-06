

% function [PB_Grid_Num, JacobianMatrix] = PB_Grid_Num(rGrid, rp1Grid, rm1Grid, rpNGrid, rmNGrid, ParaGrid, ParaSystem)
function [PB_Grid_Num] = PB_Grid_Num(rGrid, ParaGrid, ParaSystem_ini)

    global ParaSystem
    ParaSystem = ParaSystem_ini;
    
    Beta = ParaSystem.Beta;

%     ParaSystem_ini = ParaSystem;
    ThetaGrid = ParaGrid.ThetaGrid;
    PhiGrid = ParaGrid.PhiGrid;

    DeltaTheta = ParaGrid.DeltaTheta;
    DeltaPhi = ParaGrid.DeltaPhi;

    NoseDirection = ParaSystem.NoseDirection;

    Nose_r0 = norm(ParaSystem.Nose.KSMU_Rp .* ParaSystem.Rp ./ ParaSystem.r0);

    %% TEST
%     rGrid(:, 1) = Nose_r0;
    
    %% Solar Wind Pressure

%{
    
% Partial Derivatives
    rGrid_mNphi = rmNGrid;

    rGrid_pNphi = rpNGrid;
    rRight = 2 * rGrid(:, end) - rGrid_mNphi(:, end);
    rGrid_pNphi(:, end) = rRight;

    rGrid_p1 = rp1Grid;

    rGrid_m1 = rm1Grid;

    drdthetaGrid_BWD = (rGrid - rGrid_mNphi) ./ (DeltaTheta);
    drdthetaGrid_FWD = (rGrid_pNphi - rGrid) ./ (DeltaTheta);
    drdthetaGrid_Central = (rGrid_pNphi - rGrid_mNphi) ./ (2*DeltaTheta);
% 
    C = 1;
    drdthetaGrid = C*drdthetaGrid_Central + (1-C)*drdthetaGrid_FWD;
% 
%     C = 1;
%     drdthetaGrid = C*drdthetaGrid_BWD + (1-C)*drdthetaGrid_FWD;

%     drdthetaGrid = drdthetaGrid_Central;
    % drdthetaGrid(:,1) = drdthetaGrid_FWD(:,1);
    % drdthetaGrid(:,2:end) = drdthetaGrid_Central(:,2:end);
%     drdthetaGrid = drdthetaGrid_FWD;

%     drdthetaGrid(:, 1) = -rGrid(:, 1) .* tan( acos(dot(ParaSystem.NoseDirection, [1; 0; 0])) );

    drdphiGrid = (rGrid_p1 - rGrid_m1) ./ (2*DeltaPhi);
%     figure;
%     pcolor(Psw)
%     colorbar
     
    nGrid_er = 0*rGrid + 1;
    nGrid_eTheta = ( -1./rGrid) .* drdthetaGrid;
    nGrid_ePhi = ( -1./ ( rGrid.*sin(ThetaGrid) ) ) .* drdphiGrid;
    nGrid_ePhi(:, 1) = 0;
    
    nGrid_ePhi(1, :) = 0;
    nGrid_ePhi(end, :) = 0;
    

    nY_Nose = nGrid_er .* sin(ThetaGrid) .* cos(PhiGrid) +  nGrid_eTheta .* cos(ThetaGrid) .* cos(PhiGrid) +  nGrid_ePhi .* (-sin(PhiGrid));
    nZ_Nose = nGrid_er .* sin(ThetaGrid) .* sin(PhiGrid) +  nGrid_eTheta .* cos(ThetaGrid) .* sin(PhiGrid) +  nGrid_ePhi .* (cos(PhiGrid));
    nX_Nose = nGrid_er .* cos(ThetaGrid) +  nGrid_eTheta .* (-sin(ThetaGrid));

% %% TEST
% 
%     [Y, Z, X] = sph2cart(PhiGrid, pi/2-ThetaGrid, rGrid .* ParaSystem.r0/ParaSystem.Rp);
%     [nX_Nose, nY_Nose, nZ_Nose] = surfnorm(X, Y, Z);
% %     quiver3(X, Y, Z, nX_Nose, nY_Nose, nZ_Nose, 'AutoScale', 'on', 'AutoScaleFactor', 2)
%     
% %% TEST    

%}

% Conversion of nGrid from Nose_Cartesian to KSM
    exNose = NoseDirection;

%     RotationTheta = ParaSystem.Tilt.RotationTheta;
%     RotationPhi = ParaSystem.Tilt.RotationPhi;
%     M_ini = ParaSystem.M * [0; 0; -1];
%     MTilted = RotationPhi * (RotationTheta * M_ini);
    
    MTilted = ParaSystem.MTilted;

    ezNose_temp = -(MTilted - dot(exNose, MTilted)*exNose);
    ezNose = ezNose_temp ./ norm(ezNose_temp);

    eyNose = cross(ezNose, exNose);

    ParaSystem.P_CartNose2KSM  = [exNose, eyNose, ezNose];
    P_CartNose2KSM = ParaSystem.P_CartNose2KSM;

    %{
    nX_KSM = P_CartNose2KSM(1,1) .* nX_Nose + P_CartNose2KSM(1,2) .* nY_Nose + P_CartNose2KSM(1,3) .* nZ_Nose;
    nY_KSM = P_CartNose2KSM(2,1) .* nX_Nose + P_CartNose2KSM(2,2) .* nY_Nose + P_CartNose2KSM(2,3) .* nZ_Nose;
    nZ_KSM = P_CartNose2KSM(3,1) .* nX_Nose + P_CartNose2KSM(3,2) .* nY_Nose + P_CartNose2KSM(3,3) .* nZ_Nose;

    
    n_Norm = sqrt( nX_KSM.^2 + nY_KSM.^2 + nZ_KSM.^2 );
    nX_KSM = nX_KSM ./ n_Norm;
    nY_KSM = nY_KSM ./ n_Norm;
    nZ_KSM = nZ_KSM ./ n_Norm;
    
    nX_KSM(:, 1) = 1;
    nY_KSM(:, 1) = 0;
    nZ_KSM(:, 1) = 0;    
    
    nY_KSM(1, :) = 0;
    nY_KSM(end, :) = 0;
%}    
    
%% TEST

    [Y, Z, X] = sph2cart(PhiGrid, pi/2-ThetaGrid, rGrid .* ParaSystem.r0/ParaSystem.Rp);

%     [Y, Z, X] = sph2cart(PhiGrid, Elevation, rGrid .* ParaSystem.r0/ParaSystem.Rp);
    X_KSM = P_CartNose2KSM(1,1) .* X + P_CartNose2KSM(1,2) .* Y + P_CartNose2KSM(1,3) .* Z;
    Y_KSM = P_CartNose2KSM(2,1) .* X + P_CartNose2KSM(2,2) .* Y + P_CartNose2KSM(2,3) .* Z;
    Z_KSM = P_CartNose2KSM(3,1) .* X + P_CartNose2KSM(3,2) .* Y + P_CartNose2KSM(3,3) .* Z;
    [nX_KSM, nY_KSM, nZ_KSM] = surfnorm(X_KSM, Y_KSM, Z_KSM);
    
%     nX_KSM(:, 1) =  mean(nX_KSM(:, 1));
%     nY_KSM(:, 1) =  0;
%     nZ_KSM(:, 1) =  mean(nZ_KSM(:, 1));
    
      
    %{
    nX_KSM(:, 1) =  mean(nX_KSM(:, 1))
    nY_KSM(:, 1) =  mean(nY_KSM(:, 1))
    nZ_KSM(:, 1) =  mean(nZ_KSM(:, 1))
    
%     Visualisation 
    figure
    quiver3(X, Y, Z, nX_KSM, nY_KSM, nZ_KSM, 'AutoScale', 'on', 'AutoScaleFactor', 2)
    axis equal
    %}
  
% ----------------------
    nX_KSM(:, 1) = 1;
    nY_KSM(:, 1) = 0;
    nZ_KSM(:, 1) = 0;    
% ----------------------

%     
%     nY_KSM(1, :) = 0;
%     nY_KSM(end, :) = 0;
    
    
    nNorm = sqrt(nX_KSM.^2 + nY_KSM .^2 + nZ_KSM.^2);
    nX_KSM = nX_KSM ./ nNorm;
    nY_KSM = nY_KSM ./ nNorm;
    nZ_KSM = nZ_KSM ./ nNorm;
    
%% TEST    
    

% Solar Wind Dynamic Pressure, KSM
    vDotn = -nX_KSM;
    Psw = (-1/2) * vDotn;


%% Magnetic Pressure
 
%{
Coarse = 2;
ParaGrid2 = ParaGrid;
ParaGrid2.DeltaTheta = Coarse*ParaGrid.DeltaTheta;
ParaGrid2.DeltaPhi = Coarse*ParaGrid.DeltaPhi;
ParaGrid2.PhiList = ParaGrid.PhiList(1:Coarse:end);
ParaGrid2.ThetaList = ParaGrid.ThetaList(1:Coarse:end);
[GridPhi, GridTheta] = ndgrid(ParaGrid2.PhiList, ParaGrid2.ThetaList);
ParaGrid2.PhiGrid = GridPhi;
ParaGrid2.ThetaGrid = GridTheta;
rGrid2 = rGrid(1:Coarse:end, 1:Coarse:end);

BShielding = ShieldingField(rGrid2, ParaGrid2, ParaSystem);

ParaSystem.BShielding = BShielding;
%} 

% Dipole Field + Shielding Field
[BField] = MagField(rGrid, ParaGrid, ParaSystem);

% BField.X(1, :) = BField.X(end, :);
% BField.Y(1, :) = BField.Y(end, :);
% BField.Z(1, :) = BField.Z(end, :);


% Wrapped in MagField function
%{
% Magnetic Pressure: Rotated Magnetic Moment, KSM
    MGridX = 0*rGrid + MTilted(1);
    MGridY = 0*rGrid + MTilted(2);
    MGridZ = 0*rGrid + MTilted(3);


% Magnetic Pressure: Position of Nose, KSM
    XGrid_Norm = 1 .* cos(ThetaGrid);
    YGrid_Norm = 1 .* sin(ThetaGrid) .* cos(PhiGrid);
    ZGrid_Norm = 1 .* sin(ThetaGrid) .* sin(PhiGrid);
    
    % TEST %
    XGrid_Norm = P_CartNose2KSM(1,1) .* XGrid_Norm + P_CartNose2KSM(1,2) .* YGrid_Norm + P_CartNose2KSM(1,3) .* ZGrid_Norm;
    YGrid_Norm = P_CartNose2KSM(2,1) .* XGrid_Norm + P_CartNose2KSM(2,2) .* YGrid_Norm + P_CartNose2KSM(2,3) .* ZGrid_Norm;
    ZGrid_Norm = P_CartNose2KSM(3,1) .* XGrid_Norm + P_CartNose2KSM(3,2) .* YGrid_Norm + P_CartNose2KSM(3,3) .* ZGrid_Norm;
    % TEST %
    
    YGrid_Norm(1, :) = 0;
    YGrid_Norm(end, :) = 0;

    MdotEr = MTilted(1) .* XGrid_Norm + MTilted(2) .* YGrid_Norm + MTilted(3) .* ZGrid_Norm;

    BFieldX = (1./rGrid.^3) .* ( 3*MdotEr .* XGrid_Norm - MGridX );
    BFieldY = (1./rGrid.^3) .* ( 3*MdotEr .* YGrid_Norm - MGridY );
    BFieldZ = (1./rGrid.^3) .* ( 3*MdotEr .* ZGrid_Norm - MGridZ );
%}

    nCrossB_X = nY_KSM .* BField.Z - nZ_KSM .* BField.Y;
    nCrossB_Y = nZ_KSM .* BField.X - nX_KSM .* BField.Z;
    nCrossB_Z = nX_KSM .* BField.Y - nY_KSM .* BField.X;

% Magnetic Pressure
    Pmag = sqrt(nCrossB_X.^2 + nCrossB_Y.^2 + nCrossB_Z.^2) * sqrt(1+Beta);

%% Fringing Field

%{
    % Point M on the MP, at which we want to compute Bf(M), KSM/Rp
    X_M_Nose = rGrid .* cos(ThetaGrid) .* ParaSystem.r0 ./ ParaSystem.Rp;
    Y_M_Nose = rGrid .* sin(ThetaGrid) .* cos(PhiGrid) .* ParaSystem.r0 ./ ParaSystem.Rp;
    Z_M_Nose = rGrid .* sin(ThetaGrid) .* sin(PhiGrid) .* ParaSystem.r0 ./ ParaSystem.Rp;

    X_M_Half = P_CartNose2KSM(1,1) .* X_M_Nose + P_CartNose2KSM(1,2) .* Y_M_Nose + P_CartNose2KSM(1,3) .* Z_M_Nose;
    Y_M_Half = P_CartNose2KSM(2,1) .* X_M_Nose + P_CartNose2KSM(2,2) .* Y_M_Nose + P_CartNose2KSM(2,3) .* Z_M_Nose;
    Z_M_Half = P_CartNose2KSM(3,1) .* X_M_Nose + P_CartNose2KSM(3,2) .* Y_M_Nose + P_CartNose2KSM(3,3) .* Z_M_Nose;
    
    X_M = vertcat(X_M_Half, X_M_Half);
    Y_M = vertcat(Y_M_Half, -Y_M_Half);
    Z_M = vertcat(Z_M_Half, Z_M_Half);

    % Contributions of Point P on the MP: dBp(M): Current Density and dSp
    Jp_X = (-1/ParaSystem.mu0) .* 2 .* nCrossB_X;
    Jp_Y = (-1/ParaSystem.mu0) .* 2 .* nCrossB_Y;
    Jp_Z = (-1/ParaSystem.mu0) .* 2 .* nCrossB_Z;
    
    dSp_Half = ( rGrid ).^2 .* sin(ParaGrid.ThetaGrid) .* ParaGrid.DeltaTheta .* ParaGrid.DeltaPhi;
    dSp = vertcat(dSp_Half, dSp_Half);
    
    % Vector PM

    X_List = X_M(:);
    Y_List = Y_M(:);
    Z_List = Z_M(:);
    N = length(X_List);
    
    for k = 1:N
        
        XM = X_List(k);
        YM = Y_List(k);
        ZM = Z_List(k);
        
        PM_X = (0 .* X_M + XM) - X_M;
        PM_Y = (0 .* Y_M + YM) - Y_M;
        PM_Z = (0 .* Z_M + ZM) - Z_M;
        
        PM = sqrt( PM_X.^2 + PM_Y.^2 + PM_Z.^2);
        JpdSp_X = Jp_X .* dSp;
        JpdSp_Y = Jp_Y .* dSp;
        JpdSp_Z = Jp_Z .* dSp;
        
    end
    
    %}
    %%%%%%%
    %%%%%%%
    %%%%%%% -> REPRENDS ICI: PROBLEME: We need to sample the ENTIRE
    %%%%%%% surface, we only have half in the computation of density
    %%%%%%% currents j... What to do?
    
    
    %% Pressure Balance
    PressureDifference = (Psw - Pmag);
    PressureAverage = (1/2)*(Psw + Pmag);
    PB_Grid_Num = ( (PressureDifference) ./  PressureAverage );
% 
    PB_Grid_Num(1, :) = 2*PB_Grid_Num(3,:) - PB_Grid_Num(2,:);
    PB_Grid_Num(end, :) = 2*PB_Grid_Num(end-2,:) - PB_Grid_Num(end-1,:);
% % 

%     PB_Grid_Num(1, :) = 2*PB_Grid_Num(2,:) - PB_Grid_Num(3,:);
%     PB_Grid_Num(end, :) = 2*PB_Grid_Num(end-1,:) - PB_Grid_Num(end-2,:);
% % %     PB_Grid_Num(:, 1) = 0;

%     PB_Grid_Num(:, end) = 2*PB_Grid_Num(:, end-1) - PB_Grid_Num(:, end-2);




    %     PB_Grid_Num = sum(abs(PB_Grid_Num(:)));
%     PB_Grid_Num = (PressureDifference);
% figure;
% pcolor(nCrossB_Z)
% colorbar

%{
% Estimating Normal Component of Field along MP Surface

BTot = sqrt( BField.X.^2 + BField.Y.^2 + BField.Z.^2 );
BNormal = (nX_KSM .* BField.X + nY_KSM .* BField.Y + nZ_KSM .* BField.Z);
BTangent = sqrt( (nY_KSM .* BField.Z - nZ_KSM .* BField.Y) .^2 + (nZ_KSM .* BField.X - nX_KSM .* BField.Z) .^2 + (nX_KSM .* BField.Y - nY_KSM .* BField.X) .^2 );
BDiff = 100 .* (abs(BNormal)) ./ BTot;
figure;
pcolor(BDiff)
colorbar
caxis([0 100])

mean(BDiff(:))


[Y, Z, X] = sph2cart(PhiGrid, pi/2-ThetaGrid, rGrid);
figure;
   quiver3(X, Y, Z, BField.X, BField.Y, BField.Z) 
axis equal


BNorm = sqrt( BField.X.^2 + BField.Y.^2 + BField.Z.^2 );
BNormal = (nX_KSM .* BField.X + nY_KSM .* BField.Y + nZ_KSM .* BField.Z) ./ BNorm;
figure;
pcolor(BNormal)
colorbar

BNorm = sqrt( BDipole.X.^2 + BDipole.Y.^2 + BDipole.Z.^2 );
BNormal = (nX_KSM .* BDipole.X + nY_KSM .* BDipole.Y + nZ_KSM .* BDipole.Z) ./ BNorm;
figure;
pcolor(BNormal)
colorbar

figure;
    [Y, Z, X] = sph2cart(PhiGrid, pi/2 - ThetaGrid, rGrid .* ParaSystem.r0 ./ ParaSystem.Rp);
    figure;
    hold on
        PlotSurf = surf(X, Y, Z, abs(BNormal) ./ BNorm);
        axis equal
        set(PlotSurf, 'FaceAlpha', 0.3)
        colorbar
        caxis([0 2])

quiver3(X, Y, Z, BDipole.X, BDipole.Y, BDipole.Z, 'AutoScaleFactor', 3 )
quiver3(X, Y, Z, BField.X, BField.Y, BField.Z, 'AutoScaleFactor', 3 )



%}

% 
% if nargout > 1
%     epsilon = ParaGrid.epsilon;
%     JacobianMatrix = Jacobian(rGrid, epsilon, ParaGrid, ParaSystem);
% end

end