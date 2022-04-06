
function [PB_Grid_A2, JacobianMatrix] = PB_Grid_A2(r, drdtheta, ParaGrid, ParaSystem)


    ThetaGrid = ParaGrid.ThetaGrid .';
    PhiGrid = ParaGrid.PhiGrid.';

    DeltaTheta = ParaGrid.DeltaTheta;

    NoseDirection = ParaSystem.NoseDirection;

    Nose_r0 = norm(ParaSystem.Nose.KSMU_Rp .* ParaSystem.Rp ./ ParaSystem.r0);


%% Solar Wind Pressure
    
     
    nGrid_er = 0*r + 1;
    nGrid_eTheta = ( -1./r) .* drdtheta;

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

  
    

    n_Norm = sqrt( nX_KSM.^2 + nY_KSM.^2 + nZ_KSM.^2 );
    nX_KSM = nX_KSM ./ n_Norm;
    nY_KSM = nY_KSM ./ n_Norm;
    nZ_KSM = nZ_KSM ./ n_Norm;
    
    nY_KSM = 0 .* nX_KSM;
    
    nX_KSM(1) = 1;
    nY_KSM(1) = 0;
    nZ_KSM(1) = 0;
    
% Solar Wind Dynamic Pressure, KSM
    vDotn = -nX_KSM;
    Psw = (-1/2) * vDotn;


%% Magnetic Pressure

% Magnetic Pressure: Rotated Magnetic Moment, KSM
    MGridX = 0*r + MTilted(1);
    MGridY = 0*r + MTilted(2);
    MGridZ = 0*r + MTilted(3);


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

    BFieldX = (1./r.^3) .* ( 3*MdotEr .* XGrid_Norm - MGridX );
    BFieldY = (1./r.^3) .* ( 3*MdotEr .* YGrid_Norm - MGridY );
    BFieldZ = (1./r.^3) .* ( 3*MdotEr .* ZGrid_Norm - MGridZ );

    nCrossB_X = nY_KSM .* BFieldZ - nZ_KSM .* BFieldY;
    nCrossB_Y = nZ_KSM .* BFieldX - nX_KSM .* BFieldZ;
    nCrossB_Z = nX_KSM .* BFieldY - nY_KSM .* BFieldX;

% Magnetic Pressure
    Pmag = sqrt(nCrossB_X.^2 + nCrossB_Y.^2 + nCrossB_Z.^2);


%% Pressure Balance
    PressureDifference = (Psw - Pmag);
    PressureAverage = (1/2)*(Psw + Pmag);
    PB_Grid_A2 = ( (PressureDifference) ./  PressureAverage );

%     PB_Grid_Num_NMM = (PressureDifference);


if nargout > 1
    epsilon = ParaGrid.epsilon;
    JacobianMatrix = Jacobian(r, epsilon, ParaGrid, ParaSystem);
end

end