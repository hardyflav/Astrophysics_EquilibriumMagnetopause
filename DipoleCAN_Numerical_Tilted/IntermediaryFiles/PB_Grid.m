
% function [PB_Grid, JacobianMatrix] = PB_Grid(rGuess, ParaGrid, ParaSystem)

function [PB_Grid] = PB_Grid(rGuess, ParaGrid, ParaSystem_ini)

    global ParaSystem
    ParaSystem = ParaSystem_ini;
    
    ThetaGrid = ParaGrid.ThetaGrid;
    PhiGrid = ParaGrid.PhiGrid;

    DeltaTheta = ParaGrid.DeltaTheta;
    DeltaPhi = ParaGrid.DeltaPhi;

    NoseDirection = ParaSystem.NoseDirection;

    Nose_r0 = norm(ParaSystem.Nose.KSMU_Rp .* ParaSystem.Rp ./ ParaSystem.r0);
    rGrid = rGuess;

    
%% Solar Wind Pressure
%     rGrid = reshape(rGuess, [length(ParaGrid.PhiList), length(ParaGrid.ThetaList)]);
%     rGrid(:,1) = norm(Nose_r0);
%{    
% Partial Derivatives
    rGrid_mNphi = rGrid(:, 1:end-1);
    rGrid_mNphi = horzcat( 0*rGrid(:, 1) + Nose_r0, rGrid_mNphi);

    rGrid_pNphi = rGrid(:, 2:end);
%     rRight = 2 * rGrid(:, end) - rGrid_mNphi(:, end);
    rRight = 2 * rGrid(:, end) - rGrid_mNphi(:, end);
    rGrid_pNphi = horzcat( rGrid_pNphi, rRight);

%     rGrid_p1 = rGrid(2:end, :);
%     rGrid_p1 = vertcat(rGrid_p1, rGrid(2, :));
% 
%     rGrid_m1 = rGrid(1:end-1, :);
%     rGrid_m1 = vertcat( rGrid(end-1, :), rGrid_m1 );

    rGrid_p1 = rGrid(2:end, :);
    rGrid_p1 = vertcat(rGrid_p1, rGrid(end-1, :));

    rGrid_m1 = rGrid(1:end-1, :);
    rGrid_m1 = vertcat( rGrid(2, :), rGrid_m1 );

     PB_Grid = PB_Grid_Num(rGrid, rGrid_p1, rGrid_m1, rGrid_pNphi, rGrid_mNphi, ParaGrid, ParaSystem)  ;
%}
    ParaSystem . test = 1;
    [PB_Grid] = PB_Grid_Num(rGrid, ParaGrid, ParaSystem)  ;
    %{
    
    %% TEST

    drdthetaGrid_BWD = (rGrid - rGrid_mNphi) ./ (DeltaTheta);
    drdthetaGrid_FWD = (rGrid_pNphi - rGrid) ./ (DeltaTheta);
    drdthetaGrid_Central = (rGrid_pNphi - rGrid_mNphi) ./ (2*DeltaTheta);

%     drdthetaGrid = (1/2) .* (drdthetaGrid_BWD + drdthetaGrid_FWD);
%     drdthetaGrid = drdthetaGrid_Central;
    % drdthetaGrid(:,1) = drdthetaGrid_FWD(:,1);
    % drdthetaGrid(:,2:end) = drdthetaGrid_Central(:,2:end);
    drdthetaGrid = drdthetaGrid_FWD;

    drdphiGrid = (rGrid_p1 - rGrid_m1) ./ (2*DeltaPhi);
%     figure;
%     pcolor(Psw)
%     colorbar
     
    nGrid_er = 0*rGrid + 1;
    nGrid_eTheta = (-1./rGrid) .* drdthetaGrid;
    nGrid_ePhi = ( -1./ ( rGrid.*sin(ThetaGrid) ) ) .* drdphiGrid;
    nGrid_ePhi(:, 1) = 0;

    nY_Nose = nGrid_er .* sin(ThetaGrid) .* cos(PhiGrid) +  nGrid_eTheta .* cos(ThetaGrid) .* cos(PhiGrid) +  nGrid_ePhi .* (-sin(PhiGrid));
    nZ_Nose = nGrid_er .* sin(ThetaGrid) .* sin(PhiGrid) +  nGrid_eTheta .* cos(ThetaGrid) .* sin(PhiGrid) +  nGrid_ePhi .* (cos(PhiGrid));
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

    nX_KSM(:, 1) = 1;
    nY_KSM(:, 1) = 0;
    nZ_KSM(:, 1) = 0;
    
    nX_KSM = nX_KSM ./ norm(nX_KSM);
    nY_KSM = nY_KSM ./ norm(nY_KSM);
    nZ_KSM = nZ_KSM ./ norm(nZ_KSM);

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
    

    MdotEr = MTilted(1) .* XGrid_Norm + MTilted(2) .* YGrid_Norm + MTilted(3) .* ZGrid_Norm;

    BFieldX = (1./rGrid.^3) .* ( 3*MdotEr .* XGrid_Norm - MGridX );
    BFieldY = (1./rGrid.^3) .* ( 3*MdotEr .* YGrid_Norm - MGridY );
    BFieldZ = (1./rGrid.^3) .* ( 3*MdotEr .* ZGrid_Norm - MGridZ );

    nCrossB_X = nY_KSM .* BFieldZ - nZ_KSM .* BFieldY;
    nCrossB_Y = nZ_KSM .* BFieldX - nX_KSM .* BFieldZ;
    nCrossB_Z = nX_KSM .* BFieldY - nY_KSM .* BFieldX;

% Magnetic Pressure
    Pmag = sqrt(nCrossB_X.^2 + nCrossB_Y.^2 + nCrossB_Z.^2);


%% Pressure Balance
    PressureDifference = Psw - Pmag;
    PressureAverage = (1/2)*(Psw + Pmag);
    PB_Grid = ( (PressureDifference) ./  PressureAverage );

% PB_Grid = (PressureDifference);


%}

% if nargout > 1
%     epsilon = ParaGrid.epsilon;
%     JacobianMatrix = Jacobian(rGrid, epsilon, ParaGrid, ParaSystem);
% end

end