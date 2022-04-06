
function ParaSystem_InterpPswBetaPhi = System_InterpPswBetaPhi(ParaSystem)

global DegToRad


%--- Importing CHOSEN parameters; PHI_Tilt, Psw, Rmp-dependent CAN-Disk
% Parameters
    ParaSystem_InterpPswBetaPhi.Tilt.Phi = ParaSystem.Inclination * ParaSystem.Phi;
    ParaSystem_InterpPswBetaPhi.Pswnpa = ParaSystem.Psw;
    D = 2.5;
    ParaSystem_InterpPswBetaPhi.CAN_DiskParameters = [ ParaSystem.mu0I/ParaSystem.b0 , ParaSystem.a, ParaSystem.b, D];

    
%--- System parameters
    ParaSystem_InterpPswBetaPhi.Tilt.Theta =  ParaSystem.Inclination * 27 * DegToRad; % Dipole tilt, in rad
    
    AlphaDipole = pi/2 - acos(sin(ParaSystem_InterpPswBetaPhi.Tilt.Phi)*sin(ParaSystem_InterpPswBetaPhi.Tilt.Theta));
    ParaSystem_InterpPswBetaPhi.AlphaDipole = AlphaDipole;


    ParaSystem_InterpPswBetaPhi.Tilt.RotationPhi = [ cos(ParaSystem_InterpPswBetaPhi.Tilt.Phi) -sin(ParaSystem_InterpPswBetaPhi.Tilt.Phi) 0 ;...
                                sin(ParaSystem_InterpPswBetaPhi.Tilt.Phi)  cos(ParaSystem_InterpPswBetaPhi.Tilt.Phi) 0 ;...
                                0              0             1];
                        
    ParaSystem_InterpPswBetaPhi.Tilt.RotationTheta =  [ cos(ParaSystem_InterpPswBetaPhi.Tilt.Theta) 0  sin(ParaSystem_InterpPswBetaPhi.Tilt.Theta);...
                                    0               1  0 ;...
                                    -sin(ParaSystem_InterpPswBetaPhi.Tilt.Theta) 0  cos(ParaSystem_InterpPswBetaPhi.Tilt.Theta)];

    ParaSystem_InterpPswBetaPhi.Tilt.RotationAlphaDipole = [ 1 0                 0                 ;...
                                            0 cos(AlphaDipole) -sin(AlphaDipole)  ;...
                                            0 sin(AlphaDipole)  cos(AlphaDipole)  ];

    ParaSystem_InterpPswBetaPhi.mu0 = 4*pi*10^(-7);
    ParaSystem_InterpPswBetaPhi.Bp = 20000*10^(-9);   % Equatorial field, T
    ParaSystem_InterpPswBetaPhi.Rp = 60280*10^3;      % Planet radius, m
    ParaSystem_InterpPswBetaPhi.b0 = sqrt(2*ParaSystem_InterpPswBetaPhi.mu0*ParaSystem_InterpPswBetaPhi.Pswnpa*10^(-9));                                                      % Field scale
    ParaSystem_InterpPswBetaPhi.r0 = (2*ParaSystem_InterpPswBetaPhi.Bp*ParaSystem_InterpPswBetaPhi.Rp^3/ParaSystem_InterpPswBetaPhi.b0)^(1/3);         % Distance scale
    ParaSystem_InterpPswBetaPhi.M = ParaSystem_InterpPswBetaPhi.Bp*ParaSystem_InterpPswBetaPhi.Rp^3 / (ParaSystem_InterpPswBetaPhi.b0*ParaSystem_InterpPswBetaPhi.r0^3);
    
    ParaSystem_InterpPswBetaPhi.M_ini = ParaSystem_InterpPswBetaPhi.M * [0; 0; -1];
    ParaSystem_InterpPswBetaPhi.MTilted = ParaSystem_InterpPswBetaPhi.Tilt.RotationPhi * (ParaSystem_InterpPswBetaPhi.Tilt.RotationTheta * ParaSystem_InterpPswBetaPhi.M_ini);
        
end