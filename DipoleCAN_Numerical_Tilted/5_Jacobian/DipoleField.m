


function BDipole = DipoleField(Grid, ParaSystem)

% ----- Extracting Grid
    X = Grid.X .* ParaSystem.Rp ./ ParaSystem.r0;
    Y = Grid.Y .* ParaSystem.Rp ./ ParaSystem.r0;
    Z = Grid.Z .* ParaSystem.Rp ./ ParaSystem.r0;

% ----- Extracting Dipole Parameters
    MTilted = ParaSystem.MTilted;
    MGridX = 0*X + MTilted(1);
    MGridY = 0*Y + MTilted(2);
    MGridZ = 0*Z + MTilted(3);

% ----  Setting Vectors
    R = sqrt( X.^2 + Y.^2 + Z.^2 );
    XNorm = X ./ R;
    YNorm = Y ./ R;
    ZNorm = Z ./ R;

    MdotEr = MTilted(1) .* XNorm + MTilted(2) .* YNorm + MTilted(3) .* ZNorm  ;

% ----  Dipole Field
    BDipole.X = (1./R.^3) .* ( 3*MdotEr .* XNorm - MGridX );
    BDipole.Y = (1./R.^3) .* ( 3*MdotEr .* YNorm - MGridY );
    BDipole.Z = (1./R.^3) .* ( 3*MdotEr .* ZNorm - MGridZ );
    
%     BDipole.XInterp = griddedInterpolant( X, Y, Z, BDipole.X, 'linear' );
%     BDipole.YInterp = griddedInterpolant( X, Y, Z, BDipole.Y, 'linear' );
%     BDipole.ZInterp = griddedInterpolant( X, Y, Z, BDipole.Z, 'linear' );
        
end