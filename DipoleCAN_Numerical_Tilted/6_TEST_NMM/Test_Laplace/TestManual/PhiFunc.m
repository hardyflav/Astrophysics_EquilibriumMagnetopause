

function PhiFunc = PhiFunc(Ak, Bk, phi)

    PhiFunc = ones( size(Ak, 1), size(Ak, 2));
    kList = ( 1:size(Ak, 3) ) .';
    
    for kTheta = 1:size(Ak, 1)
        
        for kr = 1:size(Ak, 2)
            AkVec = squeeze( Ak(kTheta, kr, :) );
            BkVec = squeeze( Bk(kTheta, kr, :) );
            
            PhiVec = AkVec .* cos(kList.*phi) + BkVec .* sin(kList.*phi) ;
            PhiFunc(kTheta, kr) = sum(PhiVec);
        end
        
        
    end
        
end