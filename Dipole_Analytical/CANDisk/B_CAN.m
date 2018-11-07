
function [B_CAN_r, B_CAN_theta, B_CAN_phi] = B_CAN(r, theta, phi)

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
    frho = @(x) sinh(x.*D) .* exp(-x.*abs(z)) ./ x  ;
    fz = frho;
else
    frho = @(x) sinh(x.*z) .* exp(-x.*D) ./ x  ;
	fz = @(x) ( 1 - exp(-x.*D).*cosh(x.*z) ) ./ x  ;
end


% For rho > R1
[IntegralrhoEast,reterr_rhoEast,evals_rhEasto] = IIPBF(frho,rho,R1,1,0,1e-13,1e-8,'JJ');
BrhoEast = sign(z) * mu0*I0 * IntegralrhoEast;

[IntegralzEast,reterr_zEast,evals_zEast] = IIPBF(fz,rho,R1,0,0,1e-13,1e-8,'JJ');
BzEast = mu0*I0 * IntegralzEast;


% For rho > R2
[IntegralrhoWest,reterr_rhoWest,evals_rhoWest] = IIPBF(frho,rho,R2,1,0,1e-13,1e-8,'JJ');
BrhoWest = sign(z) * mu0 * (-I0) * IntegralrhoWest;

[IntegralzWest,reterr_zWest,evals_zWest] = IIPBF(fz,rho,R2,0,0,1e-13,1e-8,'JJ');
BzWest = mu0 * (-I0) * IntegralzWest;


% Total field in CAN89 Coordinate system
B_CAN_rho =  BrhoEast + BrhoWest; % Contains er and ey components in MB64 coordinate system
B_CAN_z = BzEast + BzWest;

% Total field in MB64 Coordinate system
B_CAN_r = ( B_CAN_rho / sqrt( 1-(sin(theta)*sin(phi))^2 ) ) / b0;
B_CAN_Y = ( B_CAN_z  - (sin(theta)*sin(phi))/sqrt( 1-(sin(theta)*sin(phi))^2 ) * B_CAN_rho ) / b0;

B_CAN_r = B_CAN_r + sin(theta)*sin(phi)*B_CAN_Y;
B_CAN_theta = cos(theta)*sin(phi)*B_CAN_Y;
B_CAN_phi = cos(phi)*B_CAN_Y;

end






%%% POUBELLE

%{
%%% For rho > R1
f1 = @(lambda) besselj(1, lambda.*rho) .* besselj(0, lambda.*R1) .* sinh(lambda.*D) .* exp(-lambda.*abs(z)) ./ lambda  ;
Brho = z/abs(z) * mu0*I0 * integral(f1, 0, Inf, 'RelTol',1e-8,'AbsTol',1e-13)  

f2 = @(lambda) besselj(0, lambda.*rho) .* besselj(0, lambda.*R1) .* sinh(lambda.*D) .* exp(-lambda.*abs(z)) ./ lambda  ;
Bz = mu0*I0 * integral(f2, 0, Inf, 'RelTol',1e-8,'AbsTol',1e-13)  
%}