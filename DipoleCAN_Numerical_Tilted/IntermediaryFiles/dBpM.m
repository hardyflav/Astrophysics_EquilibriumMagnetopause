

% function dBpM = dBpM(rGrid, Pvec, Mvec, ParaGrid_ini, ParaSystem_ini, BField_km1, Component)

function dBpM = dBpM(rGrid, Pvec, Mvec, ParaGrid_ini, ParaSystem_ini, BField_km1)


    global ParaSystem
    ParaSystem = ParaSystem_ini;

    ParaGrid = ParaGrid_ini;
    %% Extracting Parameters
    M.X = Mvec(1) .* ParaSystem.Rp ./ ParaSystem.r0;
    M.Y = Mvec(2) .* ParaSystem.Rp ./ ParaSystem.r0;
    M.Z = Mvec(3) .* ParaSystem.Rp ./ ParaSystem.r0;
    
    P.Phi = Pvec(:,1);
    P.Theta = Pvec(:,2);


    MTilted = ParaSystem.MTilted;
    
    exNose = ParaSystem.NoseDirection;
    ezNose_temp = -(MTilted - dot(exNose, MTilted)*exNose);
    ezNose = ezNose_temp ./ norm(ezNose_temp);
    eyNose = cross(ezNose, exNose);
    ParaSystem.P_CartNose2KSM  = [exNose, eyNose, ezNose];    
    P_CartNose2KSM = ParaSystem.P_CartNose2KSM;

    
    
%% Building Entire MP Boundary, in r0
%     rSol = vertcat(rGrid, flipud( rGrid(1:end-1, :)) );
%     ParaGrid.PhiList = vertcat(ParaGrid_ini.PhiList, pi+ParaGrid_ini.PhiList(1:end-1, :) );
    rSol = vertcat(rGrid, flipud( rGrid(1:end-1, :)) );
    ParaGrid.PhiList = linspace(-pi/2, 3*pi/2, size(rSol, 1));
%     ParaGrid.PhiList = vertcat(ParaGrid_ini.PhiList, pi+ParaGrid_ini.PhiList(1:end-1, :) );
    
    [PhiGrid, ThetaGrid] = ndgrid(ParaGrid.PhiList, ParaGrid_ini.ThetaList);
    ParaGrid.PhiGrid = PhiGrid;
    ParaGrid.ThetaGrid = ThetaGrid;
    
    % Magnetic Field at P on MP: Dipole + Bf_{k-1}
    BField.X = BField_km1.XInterp( P.Phi, P.Theta );
    BField.Y = BField_km1.YInterp( P.Phi, P.Theta );
    BField.Z = BField_km1.ZInterp( P.Phi, P.Theta );
    
       %{
       % Visualisation
       figure;
       quiver3(XSurf_KSM, YSurf_KSM, ZSurf_KSM, reshape(BField.X, size(XSurf_KSM)), reshape(BField.Y, size(XSurf_KSM)), reshape(BField.Z, size(XSurf_KSM)));
       axis equal
        %}
    
    
%% Find normals at P point on Boundary for Bp(M), in r0
    [YSurf, ZSurf, XSurf] = sph2cart(ParaGrid.PhiGrid, pi/2-ParaGrid.ThetaGrid, rSol);
    XSurf_KSM = P_CartNose2KSM(1,1) .* XSurf + P_CartNose2KSM(1,2) .* YSurf + P_CartNose2KSM(1,3) .* ZSurf;
    YSurf_KSM = P_CartNose2KSM(2,1) .* XSurf + P_CartNose2KSM(2,2) .* YSurf + P_CartNose2KSM(2,3) .* ZSurf;
    ZSurf_KSM = P_CartNose2KSM(3,1) .* XSurf + P_CartNose2KSM(3,2) .* YSurf + P_CartNose2KSM(3,3) .* ZSurf;
    [nSurf_X_KSM_Val, nSurf_Y_KSM_Val, nSurf_Z_KSM_Val] = surfnorm(XSurf_KSM, YSurf_KSM, ZSurf_KSM);
%     nSurf_X_KSM_Val(:, 1) = 1;
%     nSurf_Y_KSM_Val(:, 1) = 0;
%     nSurf_Z_KSM_Val(:, 1) = 0;
    
    %{
    % TEST: linear expansion at the nose + at the NMM Joint
    nSurf_X_KSM_Val(:, 1) = 2*nSurf_X_KSM_Val(:, 2) - nSurf_X_KSM_Val(:, 3);
    nSurf_Y_KSM_Val(:, 1) = 2*nSurf_Y_KSM_Val(:, 2) - nSurf_Y_KSM_Val(:, 3);
    nSurf_Z_KSM_Val(:, 1) = 2*nSurf_Z_KSM_Val(:, 2) - nSurf_Z_KSM_Val(:, 3);
    
    nSurf_X_KSM_Val(1, :) = 2*nSurf_X_KSM_Val(2, :) - nSurf_X_KSM_Val(3, :) ;
    nSurf_Y_KSM_Val(1, :) = 2*nSurf_Y_KSM_Val(2, :) - nSurf_Y_KSM_Val(3, :) ;
    nSurf_Z_KSM_Val(1, :) = 2*nSurf_Z_KSM_Val(2, :) - nSurf_Z_KSM_Val(3, :) ;
    
    nSurf_X_KSM_Val(end, :) = 2*nSurf_X_KSM_Val(end-1, :) - nSurf_X_KSM_Val(end-2, :) ;
    nSurf_Y_KSM_Val(end, :) = 2*nSurf_Y_KSM_Val(end-1, :) - nSurf_Y_KSM_Val(end-2, :) ;
    nSurf_Z_KSM_Val(end, :) = 2*nSurf_Z_KSM_Val(end-1, :) - nSurf_Z_KSM_Val(end-2, :) ;
    
    nSurf_X_KSM_Val(1, :) = (1/2) .* ( nSurf_X_KSM_Val(1, :) + nSurf_X_KSM_Val(end, :) );
    nSurf_Y_KSM_Val(1, :) = (1/2) .* ( nSurf_Y_KSM_Val(1, :) + nSurf_Y_KSM_Val(end, :) );
    nSurf_Z_KSM_Val(1, :) = (1/2) .* ( nSurf_Z_KSM_Val(1, :) + nSurf_Z_KSM_Val(end, :) );
    
    nSurf_X_KSM_Val(end, :) = nSurf_X_KSM_Val(1, :) ;
    nSurf_Y_KSM_Val(end, :) = nSurf_Y_KSM_Val(1, :) ;
    nSurf_Z_KSM_Val(end, :) = nSurf_Z_KSM_Val(1, :) ;
       % TEST: linear expansion at the nose
%}
    
       %{
       % Visualisation
       figure;
       quiver3(XSurf_KSM, YSurf_KSM, ZSurf_KSM, nSurf_X_KSM_Val, nSurf_Y_KSM_Val, nSurf_Z_KSM_Val);
       axis equal
        %}
    
    Method = 'linear';
    
%     griddedInterpolant(ParaGrid.PhiGrid) nSurf_X_KSM_Val(ParaGrid.PhiGrid==pi/2)
    
    nP_XInterp = griddedInterpolant(ParaGrid.PhiGrid, ParaGrid.ThetaGrid, nSurf_X_KSM_Val, Method);
    nP_YInterp = griddedInterpolant(ParaGrid.PhiGrid, ParaGrid.ThetaGrid, nSurf_Y_KSM_Val, Method);
    nP_ZInterp = griddedInterpolant(ParaGrid.PhiGrid, ParaGrid.ThetaGrid, nSurf_Z_KSM_Val, Method);

    nP_X_KSM = nP_XInterp(P.Phi, P.Theta);
    nP_Y_KSM = nP_YInterp(P.Phi, P.Theta);
    nP_Z_KSM = nP_ZInterp(P.Phi, P.Theta);
        
    nP_X_KSM = nP_X_KSM ./ sqrt( nP_X_KSM.^2 + nP_Y_KSM.^2 + nP_Z_KSM.^2 );
    nP_Y_KSM = nP_Y_KSM ./ sqrt( nP_X_KSM.^2 + nP_Y_KSM.^2 + nP_Z_KSM.^2 );
    nP_Z_KSM = nP_Z_KSM ./ sqrt( nP_X_KSM.^2 + nP_Y_KSM.^2 + nP_Z_KSM.^2 );
    
    nP_X_KSM(P.Theta==0) = 0;
    nP_Y_KSM(P.Theta==0) = 0;
    nP_Z_KSM(P.Theta==0) = 0;
    
%% Find Coordinates of P point
    rInterp = griddedInterpolant(ParaGrid.PhiGrid, ParaGrid.ThetaGrid, rSol);
    rP = rInterp(P.Phi, P.Theta);
    [YP_Nose, ZP_Nose, XP_Nose] = sph2cart(P.Phi, pi/2-P.Theta, rP);
    
    XP = P_CartNose2KSM(1,1) .* XP_Nose +  P_CartNose2KSM(1,2) .* YP_Nose + P_CartNose2KSM(1,3) .* ZP_Nose ;
    YP = P_CartNose2KSM(2,1) .* XP_Nose +  P_CartNose2KSM(2,2) .* YP_Nose + P_CartNose2KSM(2,3) .* ZP_Nose ;
    ZP = P_CartNose2KSM(3,1) .* XP_Nose +  P_CartNose2KSM(3,2) .* YP_Nose + P_CartNose2KSM(3,3) .* ZP_Nose ;
    
       %{
       % Visualisation
       figure;
       scatter3(XP, YP, ZP);
       axis equal
        %}
    
%% Around each P: local contribution
    BFieldX = BField.X;
    BFieldY = BField.Y;
    BFieldZ = BField.Z;

    nCrossB_X = nP_Y_KSM .* BFieldZ - nP_Z_KSM .* BFieldY;
    nCrossB_Y = nP_Z_KSM .* BFieldX - nP_X_KSM .* BFieldZ;
    nCrossB_Z = nP_X_KSM .* BFieldY - nP_Y_KSM .* BFieldX;

    % Point P on MP, contributing to dBp(M): Current Density and dSp
    Jp_X = (-1/1) .* 2 .* nCrossB_X;
    Jp_Y = (-1/1) .* 2 .* nCrossB_Y;
    Jp_Z = (-1/1) .* 2 .* nCrossB_Z;

%     dSp = rP.^2 .* sin(P.Theta) .* ParaGrid.DeltaTheta .* ParaGrid.DeltaPhi ;
    
    %% TEST: SURFACE INTEGRAL
%     rkpN = horzcat( rSol(:, 2:end), 2*rSol(:, end) - rSol(:, end-1) );
%     rkp1 = vertcat( rSol(2:end, :),  rSol(2, :) );
%     rkm1 = vertcat( rSol(end-1, :), rSol(1:end-1, :) );
%     drdtheta = (rkpN - rSol) ./ ParaGrid.DeltaTheta;
%     drdphi = (rkp1 - rkm1) ./ ( 2.* ParaGrid.DeltaPhi );
%     
%     dSp = sqrt( drdtheta(:).^2 + drdphi(:).^2 + 1 ) .* ParaGrid.DeltaTheta .* ParaGrid.DeltaPhi ;    

%     drdtheta(:,1) = 0 ;
%     drdphi(:,1) = 0 ;

    rkpN = horzcat( rSol(:, 2:end), 2*rSol(:, end) - rSol(:, end-1) );
    rkp1 = vertcat( rSol(2:end, :),  rSol(end-1, :) );
    rkm1 = vertcat( rSol(2, :), rSol(1:end-1, :) );
    drdtheta = (rkpN - rSol) ./ ParaGrid.DeltaTheta;
    drdphi = (rkp1 - rkm1) ./ ( 2.* ParaGrid.DeltaPhi );
   
    dSp = sqrt( drdtheta(:).^2 + drdphi(:).^2 + 1 ) .* ParaGrid.DeltaTheta .* ParaGrid.DeltaPhi ;
    
    % TEST: n Conversion for partial derivatives
    P_KSM2CartNose = P_CartNose2KSM^(-1);
    nX_Nose = P_KSM2CartNose(1,1) .* nP_X_KSM + P_KSM2CartNose(1,2) .* nP_Y_KSM + P_KSM2CartNose(1,3) .* nP_Z_KSM;
    nY_Nose = P_KSM2CartNose(2,1) .* nP_X_KSM + P_KSM2CartNose(2,2) .* nP_Y_KSM + P_KSM2CartNose(2,3) .* nP_Z_KSM;
    nZ_Nose = P_KSM2CartNose(3,1) .* nP_X_KSM + P_KSM2CartNose(3,2) .* nP_Y_KSM + P_KSM2CartNose(3,3) .* nP_Z_KSM;
    
    nTheta = cos(ParaGrid.ThetaGrid(:)).*cos(ParaGrid.PhiGrid(:)) .* nY_Nose + cos(ParaGrid.ThetaGrid(:)).*sin(ParaGrid.PhiGrid(:)) .* nZ_Nose - sin(ParaGrid.ThetaGrid(:)) .* nX_Nose;
    nPhi = -sin(ParaGrid.PhiGrid(:)) .* nY_Nose + cos(ParaGrid.PhiGrid(:)) .* nZ_Nose ;
    
    drdtheta = nTheta .* rSol(:) ;
    drdphi = nPhi .* rSol(:) .* sin(ParaGrid.ThetaGrid(:));
%     dSp = sqrt( drdtheta(:).^2 + drdphi(:).^2 + 1 ) .* ParaGrid.DeltaTheta .* ParaGrid.DeltaPhi ;

    dOMdtheta_X = drdtheta(:);
    dOMdtheta_Y = rSol(:);
    dOMdtheta_Z = 0.*rSol(:) ;
    
    dOMdphi_X = drdphi(:);
    dOMdphi_Y = 0.*rSol(:);
    dOMdphi_Z = rSol(:).*sin(ParaGrid.ThetaGrid(:));
    
    CrossNorm2 = (dOMdtheta_Y .* dOMdphi_Z - dOMdtheta_Z .* dOMdphi_Y).^2 + (dOMdtheta_Z.* dOMdphi_X - dOMdtheta_X .* dOMdphi_Z).^2 + (dOMdtheta_X.*dOMdphi_Y - dOMdtheta_Y.*dOMdphi_X).^2;
    
    dSp = sqrt( CrossNorm2 ) .* ParaGrid.DeltaTheta .* ParaGrid.DeltaPhi ;
    % TEST: n Conversion for partial derivatives

    %% TEST: SURFACE INTEGRAL

    % Infinitesimal contribution from P at M: dBp(M)
    JpdSp_X = Jp_X .* dSp;
    JpdSp_Y = Jp_Y .* dSp;
    JpdSp_Z = Jp_Z .* dSp;
    
       %{
       % Visualisation
       JNorm = sqrt(JpdSp_X.^2 + JpdSp_Y.^2 + JpdSp_Z.^2)
       figure;
       hold on
       C = 1;
    
       X = XP(1:C:end) .* ParaSystem.r0./ ParaSystem.Rp;
       Y = YP(1:C:end).* ParaSystem.r0./ ParaSystem.Rp;
       Z = ZP(1:C:end).* ParaSystem.r0./ ParaSystem.Rp;
       Plot = quiver3( X, Y, Z, JpdSp_X(1:C:end), JpdSp_Y(1:C:end), JpdSp_Z(1:C:end));
       set(Plot, 'LineWidth', 2, 'Color', 'blue', 'MaxHeadSize', 5, 'AutoScaleFactor', 1)
       axis equal
    
       axis([-20 35 -46 46 -42 40])
       grid on
    

        Ylabel = ylabel('Y (R_p)');
        Xlabel = xlabel('X (R_p)');
        Zlabel = zlabel('Z (R_p)');
        ColourBackground = [0, 0, 0, 0.5];
        set(gca,'Color',ColourBackground)
        set(gcf,'color','w');
    
    
    
        Surface = surf( XSurf_KSM.* ParaSystem.r0./ ParaSystem.Rp, YSurf_KSM.* ParaSystem.r0./ ParaSystem.Rp, ZSurf_KSM.* ParaSystem.r0./ ParaSystem.Rp )
    
    
        set(Surface, 'FaceColor', 'blue', 'FaceAlpha', 0.5, 'EdgeColor', 'none');
    
       line([0 30], [0 0], [0 0]);
       line([0 0], [-30 30], [0 0]);
       line([0 0], [0 0], [-30 30]);
    
        figure
        pcolor( reshape(JpdSp_Y, [36, 10]) )
        colorbar
        %}
    
%% Propagating at M
    PM_X = M.X - XP;
    PM_Y = M.Y - YP;
    PM_Z = M.Z - ZP;
    PM = sqrt( PM_X.^2 + PM_Y.^2 + PM_Z.^2 ) ;

    
%% Building contribution dBp(M)
    dBpM_X = 0.* PM;
    dBpM_Y = 0.* PM;
    dBpM_Z = 0.* PM;

    PM_Cr = 60 .* ParaSystem.Rp ./ ParaSystem.r0;
    Cond = ( PM > 0 & PM < 10 );
    dBpM_X(Cond) = (1/(4*pi)) .* ( JpdSp_Y(Cond) .* PM_Z(Cond) - JpdSp_Z(Cond) .* PM_Y(Cond) ) ./ PM(Cond).^3;
    dBpM_Y(Cond) = (1/(4*pi)) .* ( JpdSp_Z(Cond) .* PM_X(Cond) - JpdSp_X(Cond) .* PM_Z(Cond) ) ./ PM(Cond).^3;
    dBpM_Z(Cond) = (1/(4*pi)) .* ( JpdSp_X(Cond) .* PM_Y(Cond) - JpdSp_Y(Cond) .* PM_X(Cond) ) ./ PM(Cond).^3;
%     
    Cond2 = not( Cond );
    dBpM_X(Cond2) = 0;
    dBpM_Y(Cond2) = 0;
    dBpM_Z(Cond2) = 0;
    
    dBpM = [dBpM_X, dBpM_Y, dBpM_Z];
%     if Component == 'X'
%         dBpM = dBpM_X;
%     elseif Component == 'Y'
%         dBpM = dBpM_Y;
%     else
%         dBpM = dBpM_Z;
%     end
% %         
       %{
       % Visualisation
       figure;
       quiver3(XP, YP, ZP, dBpM_X, dBpM_Y, dBpM_Z);
       axis equal
        %}

       %{
       % Visualisation
        dBpM_Z_Grid = reshape(sqrt(dBpM_X.^2 + dBpM_Y.^2 + dBpM_Z.^2 ), size(rSol));
        dBpM_Z_Grid = reshape(dBpM(:,3),  size(rSol));
        figure;
        pcolor(dBpM_Z_Grid)
        colorbar

        PM_Grid = reshape(PM, [36, 10]);
        %}

    
    end

 