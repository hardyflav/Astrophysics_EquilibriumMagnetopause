
function f = fCoeff(location, state, Interpolants, P_KSM2CartNose)

    N = 1; % Number of equations
    nr = length(location.x); % Number of columns
    f = zeros(N,nr); % Allocate f
    
    X = location.x;
    Y = location.y;
    Z = location.z;
    
    X_Nose = P_KSM2CartNose(1,1) .* X + P_KSM2CartNose(1,2) .* Y + P_KSM2CartNose(1,3) .* Z;
    Y_Nose = P_KSM2CartNose(2,1) .* X + P_KSM2CartNose(2,2) .* Y + P_KSM2CartNose(2,3) .* Z;
    Z_Nose = P_KSM2CartNose(3,1) .* X + P_KSM2CartNose(3,2) .* Y + P_KSM2CartNose(3,3) .* Z;
    
    [Phi, Elevation, r] = cart2sph(Y_Nose, Z_Nose, X_Nose);
    Theta = pi/2-Elevation;

%     SPx = ( Interpolants.B_DipCAN.XInterp(Phi, Theta) * Interpolants.nP_XInterp(Phi, Theta) );
%     SPy = ( Interpolants.B_DipCAN.YInterp(Phi, Theta) * Interpolants.nP_YInterp(Phi, Theta) );
%     SPz = ( Interpolants.B_DipCAN.ZInterp(Phi, Theta) * Interpolants.nP_ZInterp(Phi, Theta) );

    Bnx = ( Interpolants.B_DipCAN.XInterp(Phi, Theta) * location.nx );
    Bny = ( Interpolants.B_DipCAN.YInterp(Phi, Theta) * location.ny );
    Bnz = ( Interpolants.B_DipCAN.ZInterp(Phi, Theta) * location.nz );

    
    f(:) = -( Bnx + Bny + Bnz );
    
end