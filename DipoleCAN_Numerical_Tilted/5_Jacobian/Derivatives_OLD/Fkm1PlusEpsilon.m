

function Fkm1PlusEpsilon = Fkm1PlusEpsilon(rGrid, epsilon, ParaGrid, ParaSystem)

% ThetaList = (0:pi/10:pi/2);
% PhiList = (0:2*pi); 
% [PhiGrid, ThetaGrid] = ndgrid(PhiList, ThetaList);

ThetaGrid = ParaGrid.ThetaGrid;
PhiGrid = ParaGrid.PhiGrid;

DeltaTheta = ParaGrid.DeltaTheta;
DeltaPhi = ParaGrid.DeltaPhi;

NoseDirection = ParaSystem.NoseDirection;



% 
% rList = (2./(1+cos(ThetaList))).^(0.6);
% rGrid = repmat(rList, length(PhiList), 1);
% ParaSystem.Nose = 1;
% NoseDirection = ParaSystem.Nose.Direction;

% DeltaTheta = ParaGrid.DeltaTheta;
% DeltaPhi= ParaGrid.DeltaPhi;
% DeltaTheta=1*pi/180;
% DeltaPhi=1*pi/180;

Nose_r0 = norm(ParaSystem.Nose.KSMU_Rp .* ParaSystem.Rp ./ ParaSystem.r0);


%% Solar Wind Pressure
rGrid(:,1) = norm(Nose_r0);

% TEST
% % 
% PhiList = ParaGrid.PhiList;
% Position180 = find(PhiList == pi);
% Position90 = (1/2)*(Position180+1);
% Position270 = find(PhiList == 3*pi/2);
% 
% % 
% rGrid(Position90+1:Position180,:) = flipud(rGrid(1:Position90-1,:));
% rGrid(Position270+1:end,:) =  flipud(rGrid(Position180:Position270-1,:));
% rGrid(Position90, :) = (1/2) .* (rGrid(Position90-1, :) + rGrid(Position90+1, :) );
% rGrid(Position270, :) = (1/2) .* (rGrid(Position270-1, :) + rGrid(Position270+1, :) );

% % 


% TEST


% Partial Derivatives
rGrid_pNphi = rGrid(:, 2:end);
rGrid_pNphi = horzcat( rGrid_pNphi, 2*rGrid_pNphi(:, end) - rGrid_pNphi(:, end-1));

rGrid_mNphi = rGrid(:, 1:end-1);
rGrid_mNphi = horzcat( 0*rGrid_pNphi(:, 1) + Nose_r0, rGrid_mNphi);


rGrid_p1 = rGrid(2:end, :);
rGrid_p1 = vertcat(rGrid_p1, rGrid(1, :));


rGrid_m1 = rGrid(1:end-1, :);
rGrid_m1 = vertcat( rGrid(end, :), rGrid_m1 );

%% Adding Epsilon Offset
rGrid_m1 = rGrid_m1 + epsilon;
%% Adding Epsilon Offset


drdthetaGrid_BWD = (rGrid - rGrid_mNphi) ./ (DeltaTheta);
drdthetaGrid_FWD = (rGrid_pNphi - rGrid) ./ (DeltaTheta);
drdthetaGrid_Central = (rGrid_pNphi - rGrid_mNphi) ./ (2*DeltaTheta);

% drdthetaGrid = (1/2) .* (drdthetaGrid_BWD + drdthetaGrid_FWD);
% drdthetaGrid = drdthetaGrid_Central;
% drdthetaGrid(:,1) = drdthetaGrid_FWD(:,1);
% drdthetaGrid(:,2:end) = drdthetaGrid_Central(:,2:end);
drdthetaGrid = drdthetaGrid_FWD;



% TEST
% rGrid_p1 = rGrid(2:end, :);
% rGrid_p1 = vertcat(rGrid_p1, 2*rGrid(end, :) - rGrid(end-1, :) );
% 
% rGrid_m1 = rGrid(1:end-1, :);
% rGrid_m1 = vertcat( 2*rGrid(1, :) - rGrid(2, :), rGrid_m1 );
% TEST


drdphiGrid = (rGrid_p1 - rGrid_m1) ./ (2*DeltaPhi);

nGrid_er = 0*rGrid + 1;
nGrid_eTheta = (-1./rGrid) .* drdthetaGrid;
nGrid_ePhi = ( -1./ ( rGrid.*sin(ThetaGrid) ) ) .* drdphiGrid;
nGrid_ePhi(:, 1) = 0;

% TEST: Introduce Phi symmetry
% PhiList = ParaGrid.PhiList;
% Position90 = find(PhiList == pi/2);
% Position180 = find(PhiList == pi);
% Position270 = find(PhiList == 3*pi/2);
% 
% nGrid_ePhi(Position90,:) = 0;
% nGrid_ePhi(Position270,:) = 0;

% nGrid_ePhi(Position90+1:Position180,:) = nGrid_ePhi(1:Position90-1,:);
% nGrid_ePhi(Position270+1:end,:) = nGrid_ePhi(Position180:Position270-1,:);
% 
% TEST: Introduce Phi symmetry



% Conversion of nGrid from Nose_Spherical to Nose_Cartesian
% nX_Nose = nGrid_er .* sin(ThetaGrid) .* cos(PhiGrid) +  nGrid_eTheta .* cos(ThetaGrid) .* cos(PhiGrid) +  nGrid_ePhi .* (-sin(PhiGrid));
% nY_Nose = nGrid_er .* sin(ThetaGrid) .* sin(PhiGrid) +  nGrid_eTheta .* cos(ThetaGrid) .* sin(PhiGrid) +  nGrid_ePhi .* (cos(PhiGrid));
% nZ_Nose = nGrid_er .* cos(ThetaGrid) +  nGrid_eTheta .* (-sin(ThetaGrid));

% nX_Nose = nZ_Nose;
% nY_Nose = nX_Nose;
% nZ_Nose = nY_Nose;

nY_Nose = nGrid_er .* sin(ThetaGrid) .* cos(PhiGrid) +  nGrid_eTheta .* cos(ThetaGrid) .* cos(PhiGrid) +  nGrid_ePhi .* (-sin(PhiGrid));
nZ_Nose = nGrid_er .* sin(ThetaGrid) .* sin(PhiGrid) +  nGrid_eTheta .* cos(ThetaGrid) .* sin(PhiGrid) +  nGrid_ePhi .* (cos(PhiGrid));
nX_Nose = nGrid_er .* cos(ThetaGrid) +  nGrid_eTheta .* (-sin(ThetaGrid));


% Conversion of nGrid from Nose_Cartesian to KSM
exNose = NoseDirection;

MTilted = ParaSystem.Tilt.RotationPhi*( ParaSystem.Tilt.RotationTheta* [0; 0; -1]);
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


% Solar Wind Dynamic Pressure, KSM
vDotn = -nX_KSM;
Psw = (-1/2) * vDotn;


%% Magnetic Pressure

% Magnetic Pressure: Rotated Magnetic Moment, KSM
RotationTheta = ParaSystem.Tilt.RotationTheta;
RotationPhi = ParaSystem.Tilt.RotationPhi;
M_ini = ParaSystem.M * [0; 0; -1];
M = RotationPhi * (RotationTheta * M_ini);

MGridX = 0*rGrid + M(1);
MGridY = 0*rGrid + M(2);
MGridZ = 0*rGrid + M(3);


% Magnetic Pressure: Position of Nose, KSM
XGrid_Norm = 1 .* cos(ThetaGrid);
YGrid_Norm = 1 .* sin(ThetaGrid) .* cos(PhiGrid);
ZGrid_Norm = 1 .* sin(ThetaGrid) .* sin(PhiGrid);

MdotEr = M(1) .* XGrid_Norm + M(2) .* YGrid_Norm + M(3) .* ZGrid_Norm;

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
Fkm1PlusEpsilon = (PressureDifference) ./  PressureAverage;

% Fkm1PlusEpsilon = (PressureDifference);


end