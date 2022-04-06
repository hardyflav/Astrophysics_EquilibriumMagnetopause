clear variables
clc;
p = genpath('/Users/flavienhardy/Documents/Git/MagnetopauseTilt');
addpath(p);

global DegToRad
global RadToDeg
DegToRad = pi/180;
RadToDeg = 180/pi;

%{
%%
% TEST
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

PhiTiltList = [0; 30; 60; 90];

[GridPhi, GridTheta, GridPhiTilt] = ndgrid(ParaGrid.PhiList, ParaGrid.ThetaList, PhiTiltList);

GridR = 0*GridPhi;
GridR(:, :, 1) = rSol_0;
GridR(:, :, 2) = rSol_30;
GridR(:, :, 3) = rSol_60;
GridR(:, :, 4) = rSol_90;

rInterp = griddedInterpolant(GridPhi, GridTheta, GridPhiTilt, GridR, 'spline');
PhiTilt = 30;
PhiDeg = num2str(PhiTilt);

[GridPhi, GridTheta, GridPhiTilt] = ndgrid(ParaGrid.PhiList, ParaGrid.ThetaList, PhiTilt);
rSol = rInterp(GridPhi, GridTheta, GridPhiTilt);
ParaSystem.Tilt.Phi = PhiTilt * pi/180;
    ParaSystem.Tilt.RotationPhi = [ cos(ParaSystem.Tilt.Phi) -sin(ParaSystem.Tilt.Phi) 0 ;...
                                sin(ParaSystem.Tilt.Phi)  cos(ParaSystem.Tilt.Phi) 0 ;...
                                0              0             1];
                        
    AlphaDipole = pi/2 - acos(sin(ParaSystem.Tilt.Phi)*sin(ParaSystem.Tilt.Theta));
    ParaSystem.AlphaDipole = AlphaDipole;
    ParaSystem.Tilt.RotationAlphaDipole = [ 1 0                 0                 ;...
                                            0 cos(AlphaDipole) -sin(AlphaDipole)  ;...
                                            0 sin(AlphaDipole)  cos(AlphaDipole)  ];


%%
%}


% Load Magnetopause Surface
PhiDeg = 'Untilted';
if strcmp(PhiDeg, 'Untilted') == 0 
    filename = horzcat('MagnetopauseTilt/ResultSurfaces/Phi', PhiDeg , '.mat');
else
    filename = horzcat('MagnetopauseTilt/ResultSurfaces/Untilted.mat');
end
load(filename, 'rSol', 'ParaSystem', 'ParaGrid');


rSol = vertcat(rSol, rSol(2:end-1, :) );
ParaGrid.PhiList = vertcat(ParaGrid.PhiList, pi+ParaGrid.PhiList(2:end-1, :) );
[PhiGrid, ThetaGrid] = ndgrid(ParaGrid.PhiList, ParaGrid.ThetaList);
ParaGrid.PhiGrid = PhiGrid;
ParaGrid.ThetaGrid = ThetaGrid;


MX_List = (0:1/2:30);
BM_X_List = 0.*MX_List;
BM_Y_List = BM_X_List;
BM_Z_List = BM_X_List;

BDipX_M_List = BM_X_List;
BDipY_M_List = BM_X_List;
BDipZ_M_List = BM_X_List;


for kMx = 1:length(MX_List)
    M.X = MX_List(kMx);

    % Point M at which we want to compute Bf, KSM/Rp
    % M.X = 10;
    M.Y = ParaSystem.Nose.KSMU_Rp(2);
    M.Z = ParaSystem.Nose.KSMU_Rp(3);
    

    % Point P on MP, contributing to dBp(M): Position
    [Yp, Zp, Xp] = sph2cart(ParaGrid.PhiGrid, pi/2-ParaGrid.ThetaGrid, rSol);


    % Point P on MP, contributing to dBp(M): Normal Vector
    rkp1 = vertcat(rSol(2:end, :), rSol(end-1, :));
    rkm1 = vertcat( rSol(2, :), rSol(1:end-1, :));

    rRight = 2*rSol(:, end) - rSol(:, end-1);
    rkpN = horzcat(rSol(:, 2:end), rRight);
    rLeft = 0.*rRight + norm(ParaSystem.Nose.KSMU_Rp) .* ParaSystem.Rp ./ ParaSystem.r0;
    rkmnN = horzcat(rLeft, rSol(:, 1:end-1));

    drdthetaGrid = (rkpN - rkmnN) ./ (2*ParaGrid.DeltaTheta);
    drdphiGrid = (rkp1 - rkm1) ./ (2*ParaGrid.DeltaPhi);

    nP_er = 0.*Xp + 1;
    nP_etheta = ( -1./rSol) .* drdthetaGrid;
    nP_ephi = ( -1./(rSol.*sin(ParaGrid.ThetaGrid))) .* drdphiGrid;

    nP_Y_Nose = nP_er .* sin(ParaGrid.ThetaGrid) .* cos(ParaGrid.PhiGrid) +  nP_etheta .* cos(ParaGrid.ThetaGrid) .* cos(ParaGrid.PhiGrid) +  nP_ephi .* (-sin(ParaGrid.PhiGrid));
    nP_Z_Nose = nP_er .* sin(ParaGrid.ThetaGrid) .* sin(ParaGrid.PhiGrid) +  nP_etheta .* cos(ParaGrid.ThetaGrid) .* sin(ParaGrid.PhiGrid) + nP_ephi .* (cos(ParaGrid.PhiGrid));
    nP_X_Nose = nP_er .* cos(ParaGrid.ThetaGrid) +  nP_etheta .* (-sin(ParaGrid.ThetaGrid));

    exNose = ParaSystem.NoseDirection;

    RotationTheta = ParaSystem.Tilt.RotationTheta;
    RotationPhi = ParaSystem.Tilt.RotationPhi;
    M_ini = ParaSystem.M * [0; 0; -1];
    MTilted = RotationPhi * (RotationTheta * M_ini);

    ezNose_temp = -(MTilted - dot(exNose, MTilted)*exNose);
    ezNose = ezNose_temp ./ norm(ezNose_temp);

    eyNose = cross(ezNose, exNose);

    P_CartNose2KSM  = [exNose, eyNose, ezNose];

    nP_X_KSM = P_CartNose2KSM(1,1) .* nP_X_Nose + P_CartNose2KSM(1,2) .* nP_Y_Nose + P_CartNose2KSM(1,3) .* nP_Z_Nose;
    nP_Y_KSM = P_CartNose2KSM(2,1) .* nP_X_Nose + P_CartNose2KSM(2,2) .* nP_Y_Nose + P_CartNose2KSM(2,3) .* nP_Z_Nose;
    nP_Z_KSM = P_CartNose2KSM(3,1) .* nP_X_Nose + P_CartNose2KSM(3,2) .* nP_Y_Nose + P_CartNose2KSM(3,3) .* nP_Z_Nose;

    n_Norm = sqrt( nP_X_Nose.^2 + nP_Y_Nose.^2 + nP_Z_Nose.^2 );
    nP_X_KSM = nP_X_KSM ./ n_Norm;
    nP_Y_KSM = nP_Y_KSM ./ n_Norm;
    nP_Z_KSM = nP_Z_KSM ./ n_Norm;

    nP_Y_KSM(1, :) = 0;
    nP_Y_KSM(end, :) = 0;


    % Point P on MP, contributing to dBp(M): Magnetic Field
    MGridX = 0*rSol + MTilted(1);
    MGridY = 0*rSol + MTilted(2);
    MGridZ = 0*rSol + MTilted(3);

    XGrid_Norm_Nose = 1 .* cos(ParaGrid.ThetaGrid);
    YGrid_Norm_Nose = 1 .* sin(ParaGrid.ThetaGrid) .* cos(ParaGrid.PhiGrid);
    ZGrid_Norm_Nose = 1 .* sin(ParaGrid.ThetaGrid) .* sin(ParaGrid.PhiGrid);

    XGrid_Norm = P_CartNose2KSM(1,1) .* XGrid_Norm_Nose + P_CartNose2KSM(1,2) .* YGrid_Norm_Nose + P_CartNose2KSM(1,3) .* ZGrid_Norm_Nose;
    YGrid_Norm = P_CartNose2KSM(2,1) .* XGrid_Norm_Nose + P_CartNose2KSM(2,2) .* YGrid_Norm_Nose + P_CartNose2KSM(2,3) .* ZGrid_Norm_Nose;
    ZGrid_Norm = P_CartNose2KSM(3,1) .* XGrid_Norm_Nose + P_CartNose2KSM(3,2) .* YGrid_Norm_Nose + P_CartNose2KSM(3,3) .* ZGrid_Norm_Nose;

    YGrid_Norm(1, :) = 0;
    YGrid_Norm(end, :) = 0;

    MdotEr = MTilted(1) .* XGrid_Norm + MTilted(2) .* YGrid_Norm + MTilted(3) .* ZGrid_Norm;

    BFieldX = (1./rSol.^3) .* ( 3*MdotEr .* XGrid_Norm - MGridX );
    BFieldY = (1./rSol.^3) .* ( 3*MdotEr .* YGrid_Norm - MGridY );
    BFieldZ = (1./rSol.^3) .* ( 3*MdotEr .* ZGrid_Norm - MGridZ );

    nCrossB_X = nP_Y_KSM .* BFieldZ - nP_Z_KSM .* BFieldY;
    nCrossB_Y = nP_Z_KSM .* BFieldX - nP_X_KSM .* BFieldZ;
    nCrossB_Z = nP_X_KSM .* BFieldY - nP_Y_KSM .* BFieldX;

    % Point P on MP, contributing to dBp(M): Current Density and dSp
    Jp_X = (-1/ParaSystem.mu0) .* 2 .* nCrossB_X;
    Jp_Y = (-1/ParaSystem.mu0) .* 2 .* nCrossB_Y;
    Jp_Z = (-1/ParaSystem.mu0) .* 2 .* nCrossB_Z;

    dSp = rSol.^2 .* sin(ParaGrid.ThetaGrid) .* ParaGrid.DeltaTheta .* ParaGrid.DeltaPhi;


    % Vector PM
    Xp_KSM = P_CartNose2KSM(1,1) .* Xp + P_CartNose2KSM(1,2) .* Yp + P_CartNose2KSM(1,3) .* Zp;
    Yp_KSM = P_CartNose2KSM(2,1) .* Xp + P_CartNose2KSM(2,2) .* Yp + P_CartNose2KSM(2,3) .* Zp;
    Zp_KSM = P_CartNose2KSM(3,1) .* Xp + P_CartNose2KSM(3,2) .* Yp + P_CartNose2KSM(3,3) .* Zp;

    PM_X = (ParaSystem.Rp/ParaSystem.r0) .* (0.*rSol + M.X) - Xp_KSM;
    PM_Y = (ParaSystem.Rp/ParaSystem.r0) .* (0.*rSol + M.Y) - Yp_KSM;
    PM_Z = (ParaSystem.Rp/ParaSystem.r0) .* (0.*rSol + M.Y) - Zp_KSM;

    PM = sqrt( PM_X.^2 + PM_Y.^2 + PM_Z.^2) ;


    % Infinitesimal contribution from P at M: dBp(M)
    JpdSp_X = Jp_X .* dSp;
    JpdSp_Y = Jp_Y .* dSp;
    JpdSp_Z = Jp_Z .* dSp;

    dBp_X = (ParaSystem.mu0 / (4*pi)) .* ( JpdSp_Y .*  PM_Z - JpdSp_Z .* PM_Y) ./ (PM.^3) ;
    dBp_Y = (ParaSystem.mu0 / (4*pi)) .* ( JpdSp_Z .*  PM_X - JpdSp_X .* PM_Z) ./ (PM.^3) ;
    dBp_Z = (ParaSystem.mu0 / (4*pi)) .* ( JpdSp_X .*  PM_Y - JpdSp_Y .* PM_X) ./ (PM.^3) ;

    % Sum for all P on MP: B(M)
    BM_X = sum(dBp_X(:));
    BM_Y = sum(dBp_Y(:));
    BM_Z = sum(dBp_Z(:));

    BM_X_List(kMx) = BM_X;
    BM_Y_List(kMx) = BM_Y;
    BM_Z_List(kMx) = BM_Z;
    
    
    % Dipole Field at Point M for comparison
    PositionM = [M.X; M.Y; M.Z] .* ParaSystem.Rp ./ ParaSystem.r0;
    rM = norm(PositionM);
    erM = PositionM ./ rM;
    BM = (1./rM.^3) .* ( 3*dot(MTilted, erM) .* erM - MTilted );
    BDipX_M_List(kMx) = BM(1);
    BDipY_M_List(kMx) = BM(2);
    BDipZ_M_List(kMx) = BM(3);

end






Nose = ParaSystem.Nose.KSMU_Rp(1);
figure;
hold on
    Px = plot(MX_List, BM_X_List);
    Py = plot(MX_List, BM_Y_List);
    Pz = plot(MX_List, BM_Z_List);
    NosePosition = line([Nose Nose], [-2 2]);

    Opacity = 0.7;
    CMx = cbrewer2('PuRd', 10);
    CMy = cbrewer2('PuBu', 10);
    CMz = cbrewer2('OrRd', 10);
    set(Px, 'linewidth', 5, 'color', horzcat(CMx(4,:), Opacity))
    set(Py, 'linewidth', 5, 'color', horzcat(CMy(4,:), Opacity))
    set(Pz, 'linewidth', 5, 'color', horzcat(CMz(4,:), Opacity))
    set(NosePosition, 'linewidth', 2, 'Linestyle', '-.', 'color', 'white')
    axis([0 30 min([min(BM_X_List), min(BM_Y_List), min(BM_Z_List)]) - 0.1 max([max(BM_X_List), max(BM_Y_List), max(BM_Z_List)]) + 0.1])
    
     set(gcf,'color','w')
     axes = gca;
    set(axes, 'FontName'   , 'Helvetica',...
        'Color', [0 0 0 0.6],...
        'XAxisLocation', 'bottom',...
        'YAxisLocation', 'left',...
        'FontSize', 16,...
        'LineWidth'   , 1 ...
        );
    Ylabel = ylabel('B_f / b_0');
    Xlabel = xlabel('X_{KSM} / R_p');
    Title = title( horzcat('Shielding Field, Phi = ', PhiDeg), 'FontWeight', 'bold' );
    grid on
        
    % Dipole field at M for comparison
    BDip_Norm = sqrt(BDipX_M_List.^2 + BDipY_M_List.^2 + BDipZ_M_List.^2);
    PDipNorm = plot(MX_List, BDip_Norm);
    set(PDipNorm, 'linestyle', ':', 'linewidth', 3, 'color', horzcat(CMx(1,:), Opacity))

	Legend = legend([Px, Py, Pz, PDipNorm], ["B_x", "B_y", "B_z", "B_{dip}"]);
    set(Legend, 'FontName', 'Helvetica', 'LineWidth', 1, 'FontSize', 16, 'Color', 'white');
    
    P_Sum_X = plot(MX_List, BDipX_M_List+BM_X_List);
    P_Sum_Y = plot(MX_List, BDipY_M_List+BM_Y_List);
    P_Sum_Z = plot(MX_List, BDipZ_M_List+BM_Z_List);
    
    set(P_Sum_X, 'linewidth', 5, 'color', horzcat(CMx(4,:), Opacity))
    set(P_Sum_Y, 'linewidth', 5, 'color', horzcat(CMy(4,:), Opacity))
    set(P_Sum_Z, 'linewidth', 5, 'color', horzcat(CMz(4,:), Opacity))
    
