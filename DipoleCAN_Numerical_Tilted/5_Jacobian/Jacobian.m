
function Jacobian = Jacobian(rGrid, epsilon, ParaGrid, ParaSystem)

%% Estimating derivative vectors
% dFkdrk = (1/(2*epsilon)) .* ( FkPlusEpsilon(rGrid, epsilon, ParaGrid, ParaSystem) - FkPlusEpsilon(rGrid, -epsilon, ParaGrid, ParaSystem) );
% dFkdrkp1 = (1/(2*epsilon)) .* ( Fkp1PlusEpsilon(rGrid, epsilon, ParaGrid, ParaSystem) - Fkp1PlusEpsilon(rGrid, -epsilon, ParaGrid, ParaSystem) );
% dFkdrkm1 = (1/(2*epsilon)) .* ( Fkm1PlusEpsilon(rGrid, epsilon, ParaGrid, ParaSystem) - Fkm1PlusEpsilon(rGrid, -epsilon, ParaGrid, ParaSystem) );
% dFkdrkpN = (1/(2*epsilon)) .* ( FkpNPlusEpsilon(rGrid, epsilon, ParaGrid, ParaSystem) - FkpNPlusEpsilon(rGrid, -epsilon, ParaGrid, ParaSystem) );
% dFkdrkmN = (1/(2*epsilon)) .* ( FkmNPlusEpsilon(rGrid, epsilon, ParaGrid, ParaSystem) - FkmNPlusEpsilon(rGrid, -epsilon, ParaGrid, ParaSystem) );

    N = length(ParaGrid.PhiList);



%% TEST --------------

% Determining Offset Grids
    Nose_r0 = norm(ParaSystem.Nose.KSMU_Rp .* ParaSystem.Rp ./ ParaSystem.r0);
    rGrid(:,1) = norm(Nose_r0);
    
    rGrid_mNphi = rGrid(:, 1:end-1);
    rGrid_mNphi = horzcat( 0*rGrid(:, 1) + Nose_r0, rGrid_mNphi);

    rGrid_pNphi = rGrid(:, 2:end);
    rRight = 2 * rGrid_pNphi(:, end) - rGrid_pNphi(:, end-1);
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
    
dFkdrk = (1/(2*epsilon)) .* ( PB_Grid_Num(rGrid + epsilon, rGrid_p1, rGrid_m1, rGrid_pNphi, rGrid_mNphi, ParaGrid, ParaSystem) - PB_Grid_Num(rGrid - epsilon, rGrid_p1, rGrid_m1, rGrid_pNphi, rGrid_mNphi, ParaGrid, ParaSystem) );
dFkdrkp1 = (1/(2*epsilon)) .* ( PB_Grid_Num(rGrid, rGrid_p1 + epsilon, rGrid_m1, rGrid_pNphi, rGrid_mNphi, ParaGrid, ParaSystem) - PB_Grid_Num(rGrid, rGrid_p1 - epsilon, rGrid_m1, rGrid_pNphi, rGrid_mNphi, ParaGrid, ParaSystem) );
dFkdrkm1 = (1/(2*epsilon)) .* ( PB_Grid_Num(rGrid, rGrid_p1, rGrid_m1 + epsilon, rGrid_pNphi, rGrid_mNphi, ParaGrid, ParaSystem) - PB_Grid_Num(rGrid, rGrid_p1, rGrid_m1 - epsilon, rGrid_pNphi, rGrid_mNphi, ParaGrid, ParaSystem) );
dFkdrkpN = (1/(2*epsilon)) .* ( PB_Grid_Num(rGrid, rGrid_p1, rGrid_m1, rGrid_pNphi + epsilon, rGrid_mNphi, ParaGrid, ParaSystem) - PB_Grid_Num(rGrid, rGrid_p1, rGrid_m1, rGrid_pNphi - epsilon, rGrid_mNphi, ParaGrid, ParaSystem) );
dFkdrkmN = (1/(2*epsilon)) .* ( PB_Grid_Num(rGrid, rGrid_p1, rGrid_m1, rGrid_pNphi, rGrid_mNphi + epsilon, ParaGrid, ParaSystem) - PB_Grid_Num(rGrid, rGrid_p1, rGrid_m1, rGrid_pNphi, rGrid_mNphi - epsilon, ParaGrid, ParaSystem) );




%% TEST --------------


    
%% Filling in Jacobian Matrix: Diagonal
    I1 = (1:length(ParaGrid.PhiList)*length(ParaGrid.ThetaList)).';
    J1 = I1;

    V1Grid = dFkdrk;
    
    % Correcting for linear expansion on right boundary
    RightBoundary = dFkdrkpN(:, end) ;
    V1Grid(:, end) = V1Grid(:, end) + 2.* RightBoundary;
    
    V1Grid(:, 1) = 0;
    V1 = V1Grid(:);
    
    %% Comparing Jacobian with Grad_fd
%     V1 - diag(grad_fd);    
%     Diag_fd = diag(grad_fd)  ;
%     ErrorJac = abs(grad_fd);
%     figure;
%     pcolor(ErrorJac)
%     colorbar
%     caxis([0 25])
    
%% Filling in Jacobian Matrix: Diagonal + 1
    I2 = I1;
    J2 = I2 + 1;
    V2 = 0*I2;

    V2 = dFkdrkp1(:);

    I2 = I2(1:end-1);
    J2 = J2(1:end-1);
    V2 = V2(1:end-1);

%% Filling in Jacobian Matrix: Diagonal - 1
    I3 = I1;
    J3 = I3 - 1;
    V3 = 0*I3;

    V3 = dFkdrkm1(:);

    I3 = I3(2:end);
    J3 = J3(2:end);
    V3 = V3(2:end);

%% Filling in Jacobian Matrix: Diagonal + N

    I4 = I1;
    J4 = I4 + N;
    V4 = 0*I4;

    V4 = dFkdrkpN(:);


    I4 = I4(1:end-N);
    J4 = J4(1:end-N);
    V4 = V4(1:end-N);
    
    
 %% Filling in Jacobian Matrix: Diagonal - N

    I5 = I1;
    J5 = I5 - N;
    V5 = 0*I5;

    V5 = dFkdrkmN(:);

    I5 = I5(N+1:end);
    J5 = J5(N+1:end);
    V5 = V5(N+1:end);
    
%     I5Grid = reshape(I1, [length(ParaGrid.PhiList),  length(ParaGrid.ThetaList)]);
%     I5Grid_RightBoundary = I5Grid(:, end);
%     I5 = I5Grid_RightBoundary(:);
    
%     J5 = I5 - N;
%     V5 = 0*I5;
%     
%     V5 = dFkdrkmN(:);
    
    % Correcting for linear expansion on right boundary
    RightBoundary = dFkdrkpN(:, end) ;
    V5Flipped = flipud(V5);
    V5Flipped(1:N) = V5Flipped(1:N) - flipud(RightBoundary);
    V5 = flipud(V5Flipped);
%     V5 = V5 - 1.* RightBoundary;
    
    

%% Building  Jacobian Matrix
    I = vertcat(I1, I2, I3, I4, I5);
    J = vertcat(J1, J2, J3, J4, J5);
    V = vertcat(V1, V2, V3, V4, V5);

Jacobian = sparse(I, J, V);



end