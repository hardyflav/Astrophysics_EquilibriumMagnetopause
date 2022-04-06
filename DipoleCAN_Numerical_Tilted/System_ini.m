
function ParaSystem = System_ini(Planet, Inclination, Conditions_PhiTilt_rMP_Psw)

global DegToRad
global RadToDeg

    if isnan(Conditions_PhiTilt_rMP_Psw(3))
        
        PhiTilt = Inclination * Conditions_PhiTilt_rMP_Psw(1) ;
        rMP = Conditions_PhiTilt_rMP_Psw(2);
        Pswnpa = 0.*30*0.0035 + 0.*30*0.0052 + 0.*35*0.0029 + 0*25 + 0.0097 + 0*PhiRmpPsw_Interp( Conditions_PhiTilt_rMP_Psw );
        
        ParaSystem.rMP = rMP;
        ParaSystem.Tilt.Phi = PhiTilt * DegToRad;
        ParaSystem.Pswnpa = Pswnpa;
        ParaSystem.a = 0.030 * rMP + 6.1; % Rs
        ParaSystem.b = 0.63 * rMP + 4.8; % Rs
        ParaSystem.D = 2.5; % Rs
        ParaSystem.mu0I = 1.12 * rMP + 26.7; % nT
        ParaSystem.CAN_DiskParameters = [-ParaSystem.mu0I ParaSystem.a ParaSystem.b ParaSystem.D];
        
    elseif isnan(Conditions_PhiTilt_rMP_Psw(2))
        
        PhiTilt = Inclination * Conditions_PhiTilt_rMP_Psw(1) ;
        Pswnpa = Conditions_PhiTilt_rMP_Psw(3);
        rMP = 21.51 + 0.*PhiRmpPsw_Interp(Conditions_PhiTilt_rMP_Psw);
        
        ParaSystem.Tilt.Phi = PhiTilt * DegToRad ;
        ParaSystem.Pswnpa =  Pswnpa;
        ParaSystem.rMP =  rMP;
        
        ParaSystem.a = 0.030 * rMP + 6.1; % Rs
        ParaSystem.b = 0.63 * rMP + 4.8; % Rs
        ParaSystem.D = 2.5; % Rs
        ParaSystem.mu0I = 1.12 * rMP + 26.7; % nT
        ParaSystem.CAN_DiskParameters = [-ParaSystem.mu0I ParaSystem.a ParaSystem.b ParaSystem.D];

    else
        disp('Error: Specify two out of three parameters.')        
    end

%--- System parameters
    ParaSystem.('Planet') = Planet;
    ParaSystem.Tilt.Theta =  Inclination * 27 * DegToRad; % Dipole tilt, in rad
%     ParaSystem.Tilt.Phi =  Inclination * 48 * DegToRad; % Dipole tilt, in rad
    
    AlphaDipole = pi/2 - acos(sin(ParaSystem.Tilt.Phi)*sin(ParaSystem.Tilt.Theta));
    ParaSystem.AlphaDipole = (AlphaDipole);


    ParaSystem.Tilt.RotationPhi = [ cos(ParaSystem.Tilt.Phi) -sin(ParaSystem.Tilt.Phi) 0 ;...
                                sin(ParaSystem.Tilt.Phi)  cos(ParaSystem.Tilt.Phi) 0 ;...
                                0              0             1];
                        
    ParaSystem.Tilt.RotationTheta =  [ cos(ParaSystem.Tilt.Theta) 0  sin(ParaSystem.Tilt.Theta);...
                                    0               1  0 ;...
                                    -sin(ParaSystem.Tilt.Theta) 0  cos(ParaSystem.Tilt.Theta)];

    ParaSystem.Tilt.RotationAlphaDipole = [ 1 0                 0                 ;...
                                            0 cos(AlphaDipole) -sin(AlphaDipole)  ;...
                                            0 sin(AlphaDipole)  cos(AlphaDipole)  ];
%     ParaSystem.Pswnpa = 0.0092;
    ParaSystem.Beta = 0;
    ParaSystem.mu0 = 4*pi*10^(-7);
    ParaSystem.Bp = 20000*10^(-9);   % Equatorial field, T
    ParaSystem.Rp = 60280*10^3;      % Planet radius, m
    ParaSystem.b0 = sqrt(2*ParaSystem.mu0*ParaSystem.Pswnpa*10^(-9));                                                      % Field scale
    ParaSystem.r0 = (2*ParaSystem.Bp*ParaSystem.Rp^3/ParaSystem.b0)^(1/3);         % Distance scale
    ParaSystem.M = ParaSystem.Bp*ParaSystem.Rp^3 / (ParaSystem.b0*ParaSystem.r0^3);
    
    ParaSystem.M_ini = ParaSystem.M * [0; 0; -1];
    ParaSystem.MTilted = ParaSystem.Tilt.RotationPhi * (ParaSystem.Tilt.RotationTheta * ParaSystem.M_ini);
        
    ParaSystem.CAN_DiskParameters(1) = ParaSystem.CAN_DiskParameters(1) * 1.e-9 ./ ParaSystem.b0;
    
    %{
%--- CAN Disk parameters
    ParaSystem.mu0I =   ( -60.4*10^(-9) *0  - 56.6*10^(-9) /ParaSystem.b0 ) ;      % (mu0 I0 = 60.4 nT)
    ParaSystem.a =  7;      % Inner radius of the current distribution (in Rp)
    ParaSystem.b =  16;   % Outer radius of the current distribution (in Rp)
    ParaSystem.D =  2.5;      % Disk Semi-thickness (in Rp)
    ParaSystem.CAN_DiskParameters = [ParaSystem.mu0I ParaSystem.a ParaSystem.b ParaSystem.D];
    %}
end