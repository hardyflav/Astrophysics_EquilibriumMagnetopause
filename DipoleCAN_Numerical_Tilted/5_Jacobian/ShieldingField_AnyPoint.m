

function ShieldingField = ShieldingField_AnyPoint(rGrid, GridM, ParaGrid, ParaSystem_ini, BField_km1)
 
    global ParaSystem
    ParaSystem = ParaSystem_ini;

    %% Extracting Parameters
    PhiGrid = ParaGrid.PhiGrid;
    ThetaGrid = ParaGrid.ThetaGrid;
    PhiList = ParaGrid.PhiList;
    ThetaList = ParaGrid.ThetaList;
    
    PhiGrid_ini = PhiGrid;
    ThetaGrid_ini = ThetaGrid;
    PhiList_ini = PhiList;
    ThetaList_ini = ThetaList;
    
    MTilted = ParaSystem.MTilted;
    
    exNose = ParaSystem.NoseDirection;
    ezNose_temp = -(MTilted - dot(exNose, MTilted)*exNose);
    ezNose = ezNose_temp ./ norm(ezNose_temp);
    eyNose = cross(ezNose, exNose);
    ParaSystem.P_CartNose2KSM  = [exNose, eyNose, ezNose];    
    P_CartNose2KSM = ParaSystem.P_CartNose2KSM;

    
%% Building Entire MP Boundary, in r0
    rSol = vertcat(rGrid, flipud( rGrid(2:end, :)) );
    ParaGrid.PhiList = vertcat(ParaGrid.PhiList, pi+ParaGrid.PhiList(2:end, :) );
    [PhiGrid, ThetaGrid] = ndgrid(ParaGrid.PhiList, ParaGrid.ThetaList);
    ParaGrid.PhiGrid = PhiGrid;
    ParaGrid.ThetaGrid = ThetaGrid;
    
    % Magnetic Field at M on MP: Dipole + Bf_{k-1}
%     BField = MagField(rSol, ParaGrid, ParaSystem);
    BField.X = BField_km1.XInterp( PhiGrid, ThetaGrid );
    BField.Y = BField_km1.YInterp( PhiGrid, ThetaGrid );
    BField.Z = BField_km1.ZInterp( PhiGrid, ThetaGrid );
    
    
%% Prepping M anywhere, for B(M), in r0

C = 1;

    MX_List = GridM.XList(1:C:end, 1:C:end) .* ParaSystem.Rp / ParaSystem.r0;
    MY_List = GridM.YList(1:C:end, 1:C:end).* ParaSystem.Rp / ParaSystem.r0;
    MZ_List = GridM.ZList(1:C:end, 1:C:end).* ParaSystem.Rp / ParaSystem.r0;

    BM_X_List = 0.* GridM.XGrid;
    BM_Y_List = BM_X_List;
    BM_Z_List = BM_X_List;

    PhiGrid2 = PhiGrid(1:C:end, 1:C:end);
    ThetaGrid2 =  ThetaGrid(1:C:end, 1:C:end);
    rSol2 = rSol(1:C:end, 1:C:end);
    
%     BM_X_List = 0.* rSol2(:);
%     BM_Y_List = BM_X_List;
%     BM_Z_List = BM_X_List;
%     
%     [YM_Nose, ZM_Nose, XM_Nose] = sph2cart(PhiGrid2(:), pi/2-ThetaGrid2(:), rSol2(:) );
%     MX_List = P_CartNose2KSM(1,1) .* XM_Nose + P_CartNose2KSM(1,2) .* YM_Nose + P_CartNose2KSM(1,3) .* ZM_Nose;
%     MY_List = P_CartNose2KSM(2,1) .* XM_Nose + P_CartNose2KSM(2,2) .* YM_Nose + P_CartNose2KSM(2,3) .* ZM_Nose;
%     MZ_List = P_CartNose2KSM(3,1) .* XM_Nose + P_CartNose2KSM(3,2) .* YM_Nose + P_CartNose2KSM(3,3) .* ZM_Nose;
%               
    
    
%% TEST: Prepping P points on Boundary for Bp(M), in r0
        [Yp, Zp, Xp] = sph2cart(ParaGrid.PhiGrid, pi/2-ParaGrid.ThetaGrid, rSol);
        Xp_KSM = P_CartNose2KSM(1,1) .* Xp + P_CartNose2KSM(1,2) .* Yp + P_CartNose2KSM(1,3) .* Zp;
        Yp_KSM = P_CartNose2KSM(2,1) .* Xp + P_CartNose2KSM(2,2) .* Yp + P_CartNose2KSM(2,3) .* Zp;
        Zp_KSM = P_CartNose2KSM(3,1) .* Xp + P_CartNose2KSM(3,2) .* Yp + P_CartNose2KSM(3,3) .* Zp;
        [nP_X_KSM, nP_Y_KSM, nP_Z_KSM] = surfnorm(Xp_KSM, Yp_KSM, Zp_KSM);
        nP_X_KSM(:, 1) = 1;
        nP_Y_KSM(:, 1) = 0;
        nP_Z_KSM(:, 1) = 0;    
        
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

        dSp = rSol.^2 .* sin(ParaGrid.ThetaGrid) .* ParaGrid.DeltaTheta .* ParaGrid.DeltaPhi;
        
        % Infinitesimal contribution from P at M: dBp(M)
        JpdSp_X = Jp_X .* dSp;
        JpdSp_Y = Jp_Y .* dSp;
        JpdSp_Z = Jp_Z .* dSp;
    
%% Looping over every point P of the MP boundary

    for kMX = 1:length( MX_List )
        
            for kMY = 1:length( MY_List )

                for kMZ = 1:length( MZ_List )
                
%         ( kM/length( BM_X_List ) ) * 100

        M.X = MX_List(kMX);
        M.Y = MY_List(kMY);
        M.Z = MZ_List(kMZ);
%{
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
%}
        
        %{
        exNose = ParaSystem.NoseDirection;
        ezNose_temp = -(MTilted - dot(exNose, MTilted)*exNose);
        ezNose = ezNose_temp ./ norm(ezNose_temp);

        eyNose = cross(ezNose, exNose);

        P_CartNose2KSM  = [exNose, eyNose, ezNose];
        %}
%{
        nP_X_KSM = P_CartNose2KSM(1,1) .* nP_X_Nose + P_CartNose2KSM(1,2) .* nP_Y_Nose + P_CartNose2KSM(1,3) .* nP_Z_Nose;
        nP_Y_KSM = P_CartNose2KSM(2,1) .* nP_X_Nose + P_CartNose2KSM(2,2) .* nP_Y_Nose + P_CartNose2KSM(2,3) .* nP_Z_Nose;
        nP_Z_KSM = P_CartNose2KSM(3,1) .* nP_X_Nose + P_CartNose2KSM(3,2) .* nP_Y_Nose + P_CartNose2KSM(3,3) .* nP_Z_Nose;

        n_Norm = sqrt( nP_X_Nose.^2 + nP_Y_Nose.^2 + nP_Z_Nose.^2 );
        nP_X_KSM = nP_X_KSM ./ n_Norm;
        nP_Y_KSM = nP_Y_KSM ./ n_Norm;
        nP_Z_KSM = nP_Z_KSM ./ n_Norm;

        nP_Y_KSM(1, :) = 0;
        nP_Y_KSM(end, :) = 0;
%}
        %{
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
       %}
        
%{
        BFieldX = BField.X;
        BFieldY = BField.Y;
        BFieldZ = BField.Z;        
        
        nCrossB_X = nP_Y_KSM .* BFieldZ - nP_Z_KSM .* BFieldY;
        nCrossB_Y = nP_Z_KSM .* BFieldX - nP_X_KSM .* BFieldZ;
        nCrossB_Z = nP_X_KSM .* BFieldY - nP_Y_KSM .* BFieldX;

        % Point P on MP, contributing to dBp(M): Current Density and dSp
        Jp_X = (-1/ParaSystem.mu0) .* 2 .* nCrossB_X;
        Jp_Y = (-1/ParaSystem.mu0) .* 2 .* nCrossB_Y;
        Jp_Z = (-1/ParaSystem.mu0) .* 2 .* nCrossB_Z;

        dSp = rSol.^2 .* sin(ParaGrid.ThetaGrid) .* ParaGrid.DeltaTheta .* ParaGrid.DeltaPhi;
%}

        % Vector PM
%         Xp_KSM = P_CartNose2KSM(1,1) .* Xp + P_CartNose2KSM(1,2) .* Yp + P_CartNose2KSM(1,3) .* Zp;
%         Yp_KSM = P_CartNose2KSM(2,1) .* Xp + P_CartNose2KSM(2,2) .* Yp + P_CartNose2KSM(2,3) .* Zp;
%         Zp_KSM = P_CartNose2KSM(3,1) .* Xp + P_CartNose2KSM(3,2) .* Yp + P_CartNose2KSM(3,3) .* Zp;

%         M already in r0?
%         PM_X = (ParaSystem.Rp/ParaSystem.r0) .* (0.*rSol + M.X) - Xp_KSM;
%         PM_Y = (ParaSystem.Rp/ParaSystem.r0) .* (0.*rSol + M.Y) - Yp_KSM;
%         PM_Z = (ParaSystem.Rp/ParaSystem.r0) .* (0.*rSol + M.Y) - Zp_KSM;

        PM_X = (0.*rSol + M.X) - Xp_KSM;
        PM_Y = (0.*rSol + M.Y) - Yp_KSM;
        PM_Z = (0.*rSol + M.Z) - Zp_KSM;

        PM = sqrt( PM_X.^2 + PM_Y.^2 + PM_Z.^2) ;
        
%         JpdSp_X(PM == 0) = 0;
%         JpdSp_Y(PM == 0) = 0;
%         JpdSp_Z(PM == 0) = 0;

%         PM_critical = 0.05;
        
        dBp_X = 0.* dSp;
        dBp_X(PM == 0) = 0;
        dBp_X(PM > 0) = (1 / (4*pi)) .* ( JpdSp_Y(PM > 0) .*  PM_Z(PM > 0) - JpdSp_Z(PM > 0) .* PM_Y(PM > 0)) ./ (PM(PM > 0).^3) ;
%         dBp_X( PM<PM_critical ) = 0;
        
        dBp_Y = 0.* dSp;
        dBp_Y(PM == 0) = 0;
        dBp_Y(PM > 0) = (1 / (4*pi)) .* ( JpdSp_Z(PM > 0) .*  PM_X(PM > 0) - JpdSp_X(PM > 0) .* PM_Z(PM > 0)) ./ (PM(PM > 0).^3) ;
%         dBp_Y( PM<PM_critical ) = 0;
        
        dBp_Z = 0.* dSp;
        dBp_Z(PM == 0) = 0;
        dBp_Z(PM > 0) = (1 / (4*pi)) .* ( JpdSp_X(PM > 0) .*  PM_Y(PM > 0) - JpdSp_Y(PM > 0) .* PM_X(PM > 0)) ./ (PM(PM > 0).^3) ;
%         dBp_Z( PM<PM_critical ) = 0;
        
%         dBp_X = (ParaSystem.mu0 / (4*pi)) .* ( JpdSp_Y .*  PM_Z - JpdSp_Z .* PM_Y) ./ (PM.^3) ;
%         dBp_Y = (ParaSystem.mu0 / (4*pi)) .* ( JpdSp_Z .*  PM_X - JpdSp_X .* PM_Z) ./ (PM.^3) ;
%         dBp_Z = (ParaSystem.mu0 / (4*pi)) .* ( JpdSp_X .*  PM_Y - JpdSp_Y .* PM_X) ./ (PM.^3) ;
        
%     [row, col] = find(ismember(BM_X_List, max(ShieldingField.Z(:))))
%     max(dBp_Y(:))
    
        % Sum for all P on MP: B(M)
        BM_X = sum(dBp_X(:));
        BM_Y = sum(dBp_Y(:));
        BM_Z = sum(dBp_Z(:));

        BM_X_List(kMX, kMY, kMZ) = BM_X;
        BM_Y_List(kMX, kMY, kMZ) = BM_Y;
        BM_Z_List(kMX, kMY, kMZ) = BM_Z;

        
                end
    
            end
    end

%% Packaging Output

    ShieldingField.XGrid = BM_X_List;
    ShieldingField.YGrid = BM_Y_List;
    ShieldingField.ZGrid = BM_Z_List;
    
    Method = 'linear';
    ShieldingField.Xinterp = griddedInterpolant( GridM.XGrid, GridM.YGrid, GridM.ZGrid, BM_X_List, Method);
    ShieldingField.Yinterp = griddedInterpolant( GridM.XGrid, GridM.YGrid, GridM.ZGrid, BM_Y_List, Method);
    ShieldingField.Zinterp = griddedInterpolant( GridM.XGrid, GridM.YGrid, GridM.ZGrid, BM_Z_List, Method);
    

%     ShieldingField.X = reshape( BM_X_List, size(rSol2) );
%     ShieldingField.Y = reshape( BM_Y_List, size(rSol2) );
%     ShieldingField.Z = reshape( BM_Z_List, size(rSol2) );
%     
%     Method = 'linear';
%     ShieldingField.XInterp = griddedInterpolant(PhiGrid2, ThetaGrid2, ShieldingField.X, Method);
%     ShieldingField.YInterp = griddedInterpolant(PhiGrid2, ThetaGrid2, ShieldingField.Y, Method);
%     ShieldingField.ZInterp = griddedInterpolant(PhiGrid2, ThetaGrid2, ShieldingField.Z, Method);
%     
%     
% %% Reverting back to singular surface
%     ParaGrid.PhiList = PhiList_ini;
%     ParaGrid.PhiGrid = PhiGrid_ini;
%     ParaGrid.ThetaGrid = ThetaGrid_ini;
%     ParaGrid.ThetaList = ThetaList_ini;
    
end

%{

max(ShieldingField.Z(:))

figure;
pcolor(abs(ShieldingField.Y))
colorbar
%}
