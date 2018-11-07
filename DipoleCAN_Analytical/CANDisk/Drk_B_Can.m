
function [Drk_B_CAN_r, Drk_B_CAN_theta, Drk_B_CAN_phi] = Drk_B_Can(r, theta, phi)

%%% INPUT : Position of M where the field is computed - MB64 coordinate system
%%% OUTPUT: er (radial) and eY (vertical) components of CAN field in MB64 coordinate system
%           Scaled to b0.


% System parameters
Bp = 20000*10^(-9);     % Equatorial field, T
Rp = 60280*10^3;        % Planet radius, m
mu0 = 4*pi*10^(-7);
Pswnpa = 0.02;
b0 = sqrt(2*mu0*Pswnpa*10^(-9));   % Field scale
r0 = (2*Bp*Rp^3/b0)^(1/3);         % Distance scale


% CAN89 Coordinate system
z = r.*sin(theta).*sin(phi);
rho = r * sqrt( 1-(sin(theta)*sin(phi))^2 );


% Disk parameters
I0 = (60.4*10^(-9))/mu0;      % (mu0 I0 = 60.4 nT)
Rs = Rp;                      % Saturn radius, m

R1 = 8*Rs;      % Inner radius of the current distribution
R2 = 15.5*Rs;   % Outer radius of the current distribution
D = 3*Rs;       % Disk Semi-thickness

if abs(z)>D 
    frho = @(x) sinh(x.*D) .* exp(-x.*abs(z))  ;
    fz = @(x) frho  ;
else
    frho = @(x) sinh(x.*z) .* exp(-x.*D)  ;
	fz = @(x) ( 1 - exp(-x.*D).*cosh(x.*z) )  ;
end


% For rho > R1
[IntegralrhoEast1,~,~] = IIPBF(frho,rho,R1,0,0,1e-13,1e-8,'JJ');
[IntegralrhoEast2,~,~] = IIPBF(frho,rho,R1,2,0,1e-13,1e-8,'JJ');

Drk_BrhoEast = sign(z) * mu0*I0 * (IntegralrhoEast1 + IntegralrhoEast2) /2;

[IntegralzEast,~,~] = IIPBF(fz,rho,R1,1,0,1e-13,1e-8,'JJ');
Drk_BzEast = mu0*I0 * IntegralzEast * (-1);


% For rho > R2
[IntegralrhoWest1,~,~] = IIPBF(frho,rho,R2,0,0,1e-13,1e-8,'JJ');
[IntegralrhoWest2,~,~] = IIPBF(frho,rho,R2,2,0,1e-13,1e-8,'JJ');

Drk_BrhoWest = sign(z) * mu0*(-I0) * (IntegralrhoWest1 + IntegralrhoWest2) /2;

[IntegralzWest,~,~] = IIPBF(fz,rho,R2,1,0,1e-13,1e-8,'JJ');
Drk_BzWest = mu0*(-I0) * IntegralzWest * (-1);



% Total field in CAN89 Coordinate system
Drk_B_CAN_rho =  Drk_BrhoEast + Drk_BrhoWest; % Contains er and ey components in MB64 coordinate system
Drk_B_CAN_z = Drk_BzEast + Drk_BzWest;

% Total field in MB64 Coordinate system
Drk_B_CAN_r = ( Drk_B_CAN_rho / sqrt( 1-(sin(theta)*sin(phi))^2 ) ) / b0;
Drk_B_CAN_Y = ( Drk_B_CAN_z  - (sin(theta)*sin(phi))/sqrt( 1-(sin(theta)*sin(phi))^2 ) * Drk_B_CAN_rho ) / b0;

Drk_B_CAN_r = Drk_B_CAN_r + sin(theta)*sin(phi)*Drk_B_CAN_Y;
Drk_B_CAN_theta = cos(theta)*sin(phi)*Drk_B_CAN_Y;
Drk_B_CAN_phi = cos(phi)*Drk_B_CAN_Y;


end

