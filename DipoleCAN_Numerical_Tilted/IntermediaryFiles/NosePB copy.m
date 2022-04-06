
function NosePB = NosePB(r, Alpha, ParaSystem)

% Setting vectors n, v
    % Alpha = angle between Sun direction and Nose
    % From dXdAlpha = 0,
    drdAlpha = -r*tan(Alpha);
    n = [1; (1/r)*drdAlpha; 0];
    n = n ./ norm(n);
    v = [-cos(Alpha); sin(Alpha); 0];

    PressureSW = (-1/2)*dot(n,v);

% Setting Bfield
    X = r*cos(Alpha);
    Z = r*sin(Alpha) * cos(ParaSystem.AlphaDipole); % AlphaDipole: relative to (Oxz)
    Y = r*sin(Alpha) * sin(ParaSystem.AlphaDipole);
    PointCoordinates = [X; Y; Z]; % in KSMU frame
    CoordinateDetails.Type = 'Cartesian_KSM'; 
    CoordinateDetails.Units = 'r0';
    Bfield_XY = BField(CoordinateDetails, PointCoordinates, ParaSystem);
    
    eR = [cos(Alpha); sin(Alpha)*sin(ParaSystem.AlphaDipole); sin(Alpha)*cos(ParaSystem.AlphaDipole)];
    eAlpha = [cos(Alpha+pi/2); sin(Alpha+pi/2)*sin(ParaSystem.AlphaDipole); sin(Alpha+pi/2)*cos(ParaSystem.AlphaDipole)];
    ePhi = cross(eR, eAlpha);
    
    Mxy2rAlpha = [eR, eAlpha, ePhi];
    Bfield_rAlpha = (Mxy2rAlpha)^(-1) * [Bfield_XY(1); Bfield_XY(2); Bfield_XY(3)];
%     Bfield_rAlpha = vertcat(Bfield_rAlpha, 0);
    
    PressureMag = norm(cross(n, Bfield_rAlpha));
    
PBDiff = PressureMag - PressureSW;
PBAverage = (1/2)*(PressureMag + PressureSW);
% NosePB = PBDiff/PBAverage;
NosePB = PBDiff;
end