
function [PB_Grid_Num_NMM, JacobianMatrix] = PB_Grid_Num_NMM(rGrid, rpNGrid, rmNGrid, ParaGrid, ParaSystem)


    ThetaGrid = ParaGrid.ThetaGrid .';
    PhiGrid = ParaGrid.PhiGrid.';

    DeltaTheta = ParaGrid.DeltaTheta;

    NoseDirection = ParaSystem.NoseDirection;
    P_CartNose2KSM = ParaSystem.P_CartNose2KSM;

    Nose_r0 = norm(ParaSystem.Nose.KSMU_Rp .* ParaSystem.Rp ./ ParaSystem.r0);


%% Solar Wind Pressure
    
% Partial Derivatives
    rGrid_mNphi = rmNGrid;

    rGrid_pNphi = rpNGrid;
%     rRight = 2 * rGrid(:, end) - rGrid_mNphi(:, end);
%         rRight = 2 * rGrid(end) - rGrid(end-1);
%     rGrid_pNphi(end) = rRight;

    drdthetaGrid_BWD = abs(rGrid - rGrid_mNphi) ./ (DeltaTheta);
    drdthetaGrid_FWD = (rGrid_pNphi - rGrid) ./ (DeltaTheta);
    drdthetaGrid_Central = (rGrid_pNphi - rGrid_mNphi) ./ (2*DeltaTheta);

%     
%     figure;
%     hold on
%     plot(ThetaGrid .* 180/pi, drdthetaGrid_Central)
%     plot(ThetaGrid .* 180/pi, drdthetaGrid_BWD)
%     plot(ThetaGrid .* 180/pi, drdthetaGrid_FWD)
%     
    Coeff = 1;
    drdthetaGrid = Coeff .* drdthetaGrid_BWD + (1-Coeff) .* drdthetaGrid_FWD;
    
% (pi/2- atan(0.385/0.908)) .*180/pi


%     Coeff = 1/2;
%     drdthetaGrid = Coeff .* drdthetaGrid_BWD + (1-Coeff) .* drdthetaGrid_Central;
%     
%     Threshold = 0.4;
%     ThetaC = max( ThetaGrid(drdthetaGrid_Central<Threshold & ThetaGrid.* 180/pi < 85) );
%     
%     drdthetaGrid(ThetaGrid<ThetaC - 10*DeltaTheta) = drdthetaGrid_Central(ThetaGrid<ThetaC - 10*DeltaTheta);
%     drdthetaGrid(ThetaGrid == ThetaC) = drdthetaGrid_BWD(ThetaGrid == ThetaC);
%     drdthetaGrid(ThetaGrid == ThetaC + DeltaTheta) = drdthetaGrid_FWD(ThetaGrid == ThetaC + DeltaTheta);
%     drdthetaGrid(ThetaGrid > ThetaC + 10*DeltaTheta) = drdthetaGrid_Central(ThetaGrid > ThetaC + 10*DeltaTheta);
    
%     drdthetaGrid = drdthetaGrid_BWD;
%     drdthetaGrid(:,1) = drdthetaGrid_FWD(:,1);
%     drdthetaGrid(:,2:end) = drdthetaGrid_Central(:,2:end);
%     drdthetaGrid = drdthetaGrid_FWD;

%     Coeff = 1;
%     drdthetaGrid = Coeff * drdthetaGrid_Central + (1-Coeff) * drdthetaGrid_FWD;
%     drdthetaGrid(1) = -rGrid(1) .* tan( acos(ParaSystem.NoseDirection(1)) );
%     rGrid(2) = rGrid(1) + drdthetaGrid(1) .* DeltaTheta;

    %     figure;
%     pcolor(Psw)
%     colorbar
     
    nGrid_er = 0*rGrid + 1;
    nGrid_eTheta = ( -1./rGrid) .* drdthetaGrid;

    nY_Nose = nGrid_er .* sin(ThetaGrid) .* cos(PhiGrid) +  nGrid_eTheta .* cos(ThetaGrid) .* cos(PhiGrid);
    nZ_Nose = nGrid_er .* sin(ThetaGrid) .* sin(PhiGrid) +  nGrid_eTheta .* cos(ThetaGrid) .* sin(PhiGrid);
    nX_Nose = nGrid_er .* cos(ThetaGrid) +  nGrid_eTheta .* (-sin(ThetaGrid));

    

% Conversion of nGrid from Nose_Cartesian to KSM
    exNose = NoseDirection;

    RotationTheta = ParaSystem.Tilt.RotationTheta;
    RotationPhi = ParaSystem.Tilt.RotationPhi;
    M_ini = ParaSystem.M * [0; 0; -1];
    MTilted = RotationPhi * (RotationTheta * M_ini);

    ezNose_temp = -(MTilted - dot(exNose, MTilted)*exNose);
    ezNose = ezNose_temp ./ norm(ezNose_temp);

    eyNose = cross(ezNose, exNose);

    P_CartNose2KSM  = [exNose, eyNose, ezNose];

    
    nX_KSM = P_CartNose2KSM(1,1) .* nX_Nose + P_CartNose2KSM(1,2) .* nY_Nose + P_CartNose2KSM(1,3) .* nZ_Nose;
    nY_KSM = P_CartNose2KSM(2,1) .* nX_Nose + P_CartNose2KSM(2,2) .* nY_Nose + P_CartNose2KSM(2,3) .* nZ_Nose;
    nZ_KSM = P_CartNose2KSM(3,1) .* nX_Nose + P_CartNose2KSM(3,2) .* nY_Nose + P_CartNose2KSM(3,3) .* nZ_Nose;

%% TEST

    [Y, Z, X] = sph2cart(PhiGrid, pi/2-ThetaGrid, rGrid .* ParaSystem.r0/ParaSystem.Rp);
    X_KSM = P_CartNose2KSM(1,1) .* X + P_CartNose2KSM(1,2) .* Y + P_CartNose2KSM(1,3) .* Z;
    Y_KSM = P_CartNose2KSM(2,1) .* X + P_CartNose2KSM(2,2) .* Y + P_CartNose2KSM(2,3) .* Z;
    Z_KSM = P_CartNose2KSM(3,1) .* X + P_CartNose2KSM(3,2) .* Y + P_CartNose2KSM(3,3) .* Z;
    [nX_KSM, nY_KSM, nZ_KSM] = surfnorm( horzcat(X_KSM, X_KSM, X_KSM), horzcat(Y_KSM - 0.1, Y_KSM, Y_KSM + 0.1), horzcat(Z_KSM, Z_KSM, Z_KSM));

    nX_KSM = (PhiGrid(1,1) / abs(PhiGrid(1,1)) ) .* nX_KSM(:, 1);
    nY_KSM = 0.* nX_KSM;
    nZ_KSM =  (PhiGrid(1,1) / abs(PhiGrid(1,1)) ) .* nZ_KSM(:, 1);
    
    
%%
    
    n_Norm = sqrt( nX_KSM.^2 + nY_KSM.^2 + nZ_KSM.^2 );
    nX_KSM = nX_KSM ./ n_Norm;
    nY_KSM = nY_KSM ./ n_Norm;
    nZ_KSM = nZ_KSM ./ n_Norm;
    
    nY_KSM = 0 .* nX_KSM;
    
    nX_KSM(1) = 1;
    nY_KSM(1) = 0;
    nZ_KSM(1) = 0;
%      

%% TEST
%{
PhiList2 = [PhiGrid(1) - ParaGrid.DeltaPhi; PhiGrid(1); PhiGrid(1) + ParaGrid.DeltaPhi];
ThetaList2 = ParaGrid.ThetaList;
[PhiGrid2, ThetaGrid2] = ndgrid(PhiList2, ThetaList2);
rGrid2 = vertcat(rGrid.', rGrid.', rGrid.');

%     [Y, Z, X] = sph2cart(PhiGrid, pi/2-ThetaGrid, rGrid .* ParaSystem.r0/ParaSystem.Rp);
    [Y, Z, X] = sph2cart(PhiGrid2, pi/2-ThetaGrid2, rGrid2 .* ParaSystem.r0/ParaSystem.Rp);

    X_KSM = P_CartNose2KSM(1,1) .* X + P_CartNose2KSM(1,2) .* Y + P_CartNose2KSM(1,3) .* Z;
    Y_KSM = P_CartNose2KSM(2,1) .* X + P_CartNose2KSM(2,2) .* Y + P_CartNose2KSM(2,3) .* Z;
    Z_KSM = P_CartNose2KSM(3,1) .* X + P_CartNose2KSM(3,2) .* Y + P_CartNose2KSM(3,3) .* Z;
    [nX_KSM, nY_KSM, nZ_KSM] = surfnorm(X_KSM, Y_KSM, Z_KSM);
    
    nX_KSM = nX_KSM(2, :);
    nY_KSM = nY_KSM(2, :);
    nZ_KSM = nZ_KSM(2, :);

    nX_KSM(:, 1) = 1;
    nY_KSM(:, 1) = 0;
    nZ_KSM(:, 1) = 0;    
    %}
    
%% TEST

% Solar Wind Dynamic Pressure, KSM
    vDotn = -nX_KSM;
    Psw = (-1/2) * vDotn;


%% Magnetic Pressure

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

    BDipX = (1./rGrid.^3) .* ( 3*MdotEr .* XGrid_Norm - MGridX );
    BDipY = (1./rGrid.^3) .* ( 3*MdotEr .* YGrid_Norm - MGridY );
    BDipZ = (1./rGrid.^3) .* ( 3*MdotEr .* ZGrid_Norm - MGridZ );
    
    BCAN = ParaSystem.BCAN;
    BCANX = BCAN.XInterp(X_KSM, Y_KSM, Z_KSM);
    BCANY = BCAN.YInterp(X_KSM, Y_KSM, Z_KSM);
    BCANZ = BCAN.ZInterp(X_KSM, Y_KSM, Z_KSM);
    
C = ParaSystem.OnOff_CAN;
    BFieldX = BDipX + C .* BCANX;
    BFieldY = BDipY + C .* BCANY;
    BFieldZ = BDipZ + C .* BCANZ;

    nCrossB_X = nY_KSM .* BFieldZ - nZ_KSM .* BFieldY;
    nCrossB_Y = nZ_KSM .* BFieldX - nX_KSM .* BFieldZ;
    nCrossB_Z = nX_KSM .* BFieldY - nY_KSM .* BFieldX;

% Magnetic Pressure
    Pmag = sqrt(nCrossB_X.^2 + nCrossB_Y.^2 + nCrossB_Z.^2);


%% Pressure Balance
    PressureDifference = (Psw - Pmag);
    PressureAverage = (1/2)*(Psw + Pmag);
    PB_Grid_Num_NMM = ( (PressureDifference) ./  PressureAverage );

%     PB_Grid_Num_NMM = (PressureDifference);


if nargout > 1
    epsilon = ParaGrid.epsilon;
    JacobianMatrix = Jacobian(rGrid, epsilon, ParaGrid, ParaSystem);
end

end