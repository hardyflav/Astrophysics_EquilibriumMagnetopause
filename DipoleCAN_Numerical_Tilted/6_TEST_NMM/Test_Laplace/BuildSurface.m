
function Surface = BuildSurface(rSol, ParaGrid, ParaSystem)

    % TEST: Interp
    rInterp = griddedInterpolant(ParaGrid.PhiGrid, ParaGrid.ThetaGrid, rSol, 'linear');
    PhiList = linspace(ParaGrid.PhiList(1), ParaGrid.PhiList(end), length(ParaGrid.PhiList)*1);
    ThetaList = linspace(ParaGrid.ThetaList(1), ParaGrid.ThetaList(end), length(ParaGrid.ThetaList)*1);
    [PhiGrid, ThetaGrid] = ndgrid(PhiList, ThetaList);
    rVal = rInterp(PhiGrid, ThetaGrid);
    
    Elevation = pi/2 - ThetaGrid;
    [XNoseMB, YNoseMB, ZNoseMB] = sph2cart(PhiGrid, Elevation, rVal .* ParaSystem.r0 ./ ParaSystem.Rp);
    YNose = XNoseMB;
    ZNose = YNoseMB;
    XNose = ZNoseMB;
    P_CartNose2KSM = ParaSystem.P_CartNose2KSM;
    X_KSM = P_CartNose2KSM(1,1) .* XNose + P_CartNose2KSM(1,2) .* YNose + P_CartNose2KSM(1,3) .* ZNose;
    Y_KSM = P_CartNose2KSM(2,1) .* XNose + P_CartNose2KSM(2,2) .* YNose + P_CartNose2KSM(2,3) .* ZNose;
    Z_KSM = P_CartNose2KSM(3,1) .* XNose + P_CartNose2KSM(3,2) .* YNose + P_CartNose2KSM(3,3) .* ZNose;

    YNose_2 = -YNose;
    X_KSM_2 = P_CartNose2KSM(1,1) .* XNose + P_CartNose2KSM(1,2) .* YNose_2 + P_CartNose2KSM(1,3) .* ZNose;
    Y_KSM_2 = P_CartNose2KSM(2,1) .* XNose + P_CartNose2KSM(2,2) .* YNose_2 + P_CartNose2KSM(2,3) .* ZNose;
    Z_KSM_2 = P_CartNose2KSM(3,1) .* XNose + P_CartNose2KSM(3,2) .* YNose_2 + P_CartNose2KSM(3,3) .* ZNose;

    X_MP = vertcat(X_KSM, flipud(X_KSM_2));
    Y_MP = vertcat(Y_KSM, flipud(Y_KSM_2));
    Z_MP = vertcat(Z_KSM, flipud(Z_KSM_2));    
    
%     Test: Try Closing Surface
    %{
    XLast = X_MP(:,end);
    YLast = Y_MP(:,end);
    ZLast = Z_MP(:,end);
    
    X_MP2 = horzcat(X_MP, XLast);
    Y_MP2 = horzcat(Y_MP, flipud(YLast));
    Z_MP2 = horzcat(Z_MP, flipud(ZLast));
    %}

%     Option: Plot Surface
    %{ 
    figure
    surf(X_MP, Y_MP, Z_MP)
    axis equal
    %}

Surface.X_MP = X_MP;
Surface.Y_MP = Y_MP;
Surface.Z_MP = Z_MP;

    

end