
function Jacobian = Jacobian(rGrid, epsilon, ParaGrid, ParaSystem)

%% Estimating derivative vectors
dFkdrk = (1/(2*epsilon)) .* ( FkPlusEpsilon(rGrid, epsilon, ParaGrid, ParaSystem) - FkPlusEpsilon(rGrid, -epsilon, ParaGrid, ParaSystem) );
dFkdrkp1 = (1/(2*epsilon)) .* ( Fkp1PlusEpsilon(rGrid, epsilon, ParaGrid, ParaSystem) - Fkp1PlusEpsilon(rGrid, -epsilon, ParaGrid, ParaSystem) );
dFkdrkm1 = (1/(2*epsilon)) .* ( Fkm1PlusEpsilon(rGrid, epsilon, ParaGrid, ParaSystem) - Fkm1PlusEpsilon(rGrid, -epsilon, ParaGrid, ParaSystem) );
dFkdrkpN = (1/(2*epsilon)) .* ( FkpNPlusEpsilon(rGrid, epsilon, ParaGrid, ParaSystem) - FkpNPlusEpsilon(rGrid, -epsilon, ParaGrid, ParaSystem) );


%% Filling in Jacobian Matrix: Diagonal
    I1 = (1:length(ParaGrid.PhiList)*length(ParaGrid.ThetaList)).';
    J1 = I1;
    V1 = 0*I1;

    V1 = dFkdrk(:);

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
    N = length(ParaGrid.PhiList);

    I4 = I1;
    J4 = I4 + N;
    V4 = 0*I4;

    V4 = dFkdrkpN(:);


    I4 = I4(1:end-N);
    J4 = J4(1:end-N);
    V4 = V4(1:end-N);

%% Filling in Jacobian Matrix: Correction of Right Boundary
    I5 = flipud(I1);
    I5 = flipud(I5(1:N));


    J5 = I5 - 1;
    V5 = 0*I5;




%% Building  Jacobian Matrix
    I = vertcat(I1, I2, I3, I4);
    J = vertcat(J1, J2, J3, J4);
    V = vertcat(V1, V2, V3, V4);

Jacobian = sparse(I, J, V);

end