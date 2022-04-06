

function PB_Nose_BetaDetermination  = PB_Nose_BetaDetermination(r, NoseTheta, Beta, ParaSystem)

% Beta = ParaSystem.Beta;

% Solar Wind Dynamic Pressure, KSM
n = [1; 0; 0];
v = [-1; 0; 0];
Psw = (-1/2) * dot(n, v);

% Magnetic Pressure: Rotated Magnetic Moment, KSM
RotationTheta = ParaSystem.Tilt.RotationTheta;
RotationPhi = ParaSystem.Tilt.RotationPhi;
M_ini = ParaSystem.M * [0; 0; -1];
M = RotationPhi * (RotationTheta * M_ini);

% Magnetic Pressure: Position of Nose, KSM
NosePhi = pi/2-ParaSystem.AlphaDipole;
X = r*cos(NoseTheta);
Y = r*sin(NoseTheta)*cos(NosePhi);
Z = r*sin(NoseTheta)*sin(NosePhi);
Position = [X; Y; Z];
er = Position ./ norm(Position);

% Magnetic Pressure
Bfield = (1/r^3) * (    3*dot(M, er)*er - M    ); % For Dipole Field ONLY

% Adding CAN DISK to Dipole Field
BDipole = (1/r^3) * (    3*dot(M, er)*er - M    );
BCAN = ParaSystem.BCAN;

BCAN_X = BCAN.XInterp( Position(1).*ParaSystem.r0 ./ ParaSystem.Rp , Position(2).*ParaSystem.r0 ./ ParaSystem.Rp, Position(3).*ParaSystem.r0 ./ ParaSystem.Rp );
BCAN_Y = BCAN.XInterp( Position(1).*ParaSystem.r0 ./ ParaSystem.Rp , Position(2).*ParaSystem.r0 ./ ParaSystem.Rp, Position(3).*ParaSystem.r0 ./ ParaSystem.Rp );
BCAN_Z = BCAN.ZInterp( Position(1).*ParaSystem.r0 ./ ParaSystem.Rp , Position(2).*ParaSystem.r0 ./ ParaSystem.Rp, Position(3).*ParaSystem.r0 ./ ParaSystem.Rp );

C = ParaSystem.OnOff_CAN;
Bfield(1)  = BDipole(1) + C .* BCAN_X;
Bfield(2)  = BDipole(2) + C .* BCAN_Y;
Bfield(3)  = BDipole(3) + C .* BCAN_Z;

Pmag = norm( cross(n, Bfield) ) * sqrt(1+abs(Beta));

% Pressure Balance
PB_Difference = Pmag - Psw;
PB_Average = (1/2) * (Pmag + Psw);
PB_Nose_BetaDetermination = PB_Difference / PB_Average;


end