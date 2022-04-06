

clear all
clc


% Building (Theta x r) Grid
ThetaVec = linspace(0.1, pi/2, 10);
rBoundary = 20 .* (2./(1+cos(ThetaVec))).^(1/6);

rGrid = zeros(length(ThetaVec), length(ThetaVec));
ThetaGrid = zeros(length(ThetaVec), length(ThetaVec));
for k = 1:length(ThetaVec)
    rVec = linspace(1, rBoundary(k), 10);
    rGrid(:,k) = rVec;
    ThetaGrid(k,:) = ThetaVec;
end


% Initialise Phi Potential
phi = pi/2;
Ak = ones( size(ThetaGrid, 1), size(ThetaGrid, 2), 11 );
Bk = ones( size(ThetaGrid, 1), size(ThetaGrid, 2), 11 );
PhiVal = PhiFunc(Ak, Bk, phi);

% Laplacian 
Laplacian = del2(PhiVal);
Jacobian = rGrid;
Rho = rGrid .* sin(ThetaGrid);

Res = 0 .* (0:10).' ;
for k = 1:11
    Res(k) = (1./Jacobian) .* Laplacian .* PhiVal  - (k-1).*PhiVal ./ Rho ;
end





