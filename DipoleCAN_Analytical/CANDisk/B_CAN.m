
function B_Disk = B_CAN(rScaled, theta, phi)

%%% INPUT : Position of M where the field is computed - MB64 coordinate
%%% system; r is scaled to dipolar r0 (in m).

%%% OUTPUT: er (radial) and eY (vertical) components of CAN field in MB64 coordinate system
%%%         Scaled to b0.



% System parameters
Bp = 20000*10^(-9);     % Equatorial field, T
Rp = 60280*10^3;        % Planet radius, m
mu0 = 4*pi*10^(-7);
Pswnpa = 0.02;
b0 = sqrt(2*mu0*Pswnpa*10^(-9));   % Field scale
r0 = (2*Bp*Rp^3/b0)^(1/3);         % Distance scale

r = rScaled*r0/Rp;

% CAN89 Coordinate system
z = r.*sin(theta).*sin(phi);
rho = r * sqrt( 1-(sin(theta)*sin(phi))^2 );


% Disk parameters
I0 = (60.4*10^(-9))/mu0;      % (mu0 I0 = 60.4 nT)

R1 = 8;      % Inner radius of the current distribution (in Rp)
R2 = 15.5;   % Outer radius of the current distribution (in Rp)
D = 3;       % Disk Semi-thickness (in Rp)


if abs(z)>D
    frho = @(x) sinh(x.*D) .* exp(-x.*abs(z)) ./ x  ;
    fz = frho;
else
    frho = @(x) sinh(x.*z) .* exp(-x.*D) ./ x  ;
	fz = @(x) ( 1 - exp(-x.*D).*cosh(x.*z) ) ./ x  ;
end


% For rho > R1
[IntegralrhoEast,~,~] = IIPBF(frho,rho,R1,1,0,1e-13,1e-8,'JJ');
BrhoEast = sign(z) * mu0*I0 * IntegralrhoEast;

[IntegralzEast,~,~] = IIPBF(fz,rho,R1,0,0,1e-13,1e-8,'JJ');
BzEast = mu0*I0 * IntegralzEast;


% For rho > R2
[IntegralrhoWest,~,~] = IIPBF(frho,rho,R2,1,0,1e-13,1e-8,'JJ');
BrhoWest = sign(z) * mu0 * (-I0) * IntegralrhoWest;

[IntegralzWest,~,~] = IIPBF(fz,rho,R2,0,0,1e-13,1e-8,'JJ');
BzWest = mu0 * (-I0) * IntegralzWest;


% Total field in CAN89 Coordinate system
B_CAN_rho =  (BrhoEast + BrhoWest) / b0; % Contains er and ey components in MB64 coordinate system
B_CAN_z = (BzEast + BzWest) / b0;


% Total field in MB64 Coordinate system

if theta == pi/2 && phi == pi/2
    
    B_CAN_r = B_CAN_z;
    B_CAN_theta = 0;
    B_CAN_phi = 0;
    
elseif phi == 0
    
    B_CAN_r = B_CAN_rho;
    B_CAN_theta = 0;
    B_CAN_phi = B_CAN_z;
    
else
    
    [X_CAN, Y_CAN, Z_CAN] = pol2cart(0, B_CAN_rho, B_CAN_z);
    
    X_MB = Y_CAN;
    Y_MB = Z_CAN;
    Z_MB = X_CAN;
    
    [B_CAN_phi, ~, B_CAN_r] = cart2sph(X_MB, Y_MB, Z_MB);

    
    B_CAN_theta = atan2(sqrt(X_MB^2+ Y_MB^2), Z_MB);

end



B_Disk = [B_CAN_r; B_CAN_theta; B_CAN_phi];

end




%{  DE COTE
% 
% function [B_CAN_r, B_CAN_theta, B_CAN_phi] = B_CAN(rScaled, theta, phi)
% 
% %%% INPUT : Position of M where the field is computed - MB64 coordinate
% %%% system; r is scaled to r0 (in m).
% 
% %%% OUTPUT: er (radial) and eY (vertical) components of CAN field in MB64 coordinate system
% %           Scaled to b0.
% 
% 
% 
% % System parameters
% Bp = 20000*10^(-9);     % Equatorial field, T
% Rp = 60280*10^3;        % Planet radius, m
% mu0 = 4*pi*10^(-7);
% Pswnpa = 0.02;
% b0 = sqrt(2*mu0*Pswnpa*10^(-9));   % Field scale
% r0 = (2*Bp*Rp^3/b0)^(1/3);         % Distance scale
% 
% r = rScaled*r0;
% 
% % CAN89 Coordinate system
% z = r.*sin(theta).*sin(phi);
% rho = r * sqrt( 1-(sin(theta)*sin(phi))^2 );
% 
% 
% % Disk parameters
% I0 = (60.4*10^(-9))/mu0;      % (mu0 I0 = 60.4 nT)
% Rs = Rp;                      % Saturn radius, m
% 
% R1 = 8*Rs;      % Inner radius of the current distribution
% R2 = 15.5*Rs;   % Outer radius of the current distribution
% D = 3*Rs;       % Disk Semi-thickness
% 
% 
% if rho<R1
%     B_CAN_r = 0;
%     B_CAN_theta = 0;
%     B_CAN_phi = 0;
% else
% 
% if abs(z)>D 
%     frho = @(x) sinh(x.*D) .* exp(-x.*abs(z)) ./ x  ;
%     fz = frho;
% else
%     frho = @(x) sinh(x.*z) .* exp(-x.*D) ./ x  ;
% 	fz = @(x) ( 1 - exp(-x.*D).*cosh(x.*z) ) ./ x  ;
% end
% 
% 
% % For rho > R1
% [IntegralrhoEast,~,~] = IIPBF(frho,rho,R1,1,0,1e-13,1e-8,'JJ');
% BrhoEast = sign(z) * mu0*I0 * IntegralrhoEast;
% 
% [IntegralzEast,~,~] = IIPBF(fz,rho,R1,0,0,1e-13,1e-8,'JJ');
% BzEast = mu0*I0 * IntegralzEast;
% 
% 
% % For rho > R2
% [IntegralrhoWest,~,~] = IIPBF(frho,rho,R2,1,0,1e-13,1e-8,'JJ');
% BrhoWest = sign(z) * mu0 * (-I0) * IntegralrhoWest;
% 
% [IntegralzWest,~,~] = IIPBF(fz,rho,R2,0,0,1e-13,1e-8,'JJ');
% BzWest = mu0 * (-I0) * IntegralzWest;
% 
% 
% 
% % Total field in CAN89 Coordinate system
% B_CAN_rho =  BrhoEast + BrhoWest; % Contains er and ey components in MB64 coordinate system
% B_CAN_z = BzEast + BzWest;
% 
% % Total field in MB64 Coordinate system
% B_CAN_r = ( B_CAN_rho / sqrt( 1-(sin(theta)*sin(phi))^2 ) ) / b0;
% B_CAN_Y = ( B_CAN_z  - (sin(theta)*sin(phi))/sqrt( 1-(sin(theta)*sin(phi))^2 ) * B_CAN_rho ) / b0;
% 
% B_CAN_r = B_CAN_r + sin(theta)*sin(phi)*B_CAN_Y;
% B_CAN_theta = cos(theta)*sin(phi)*B_CAN_Y;
% B_CAN_phi = cos(phi)*B_CAN_Y;
% 
% end
% 
% end




% B_CAN_r = ( B_CAN_rho / sqrt( 1-(sin(theta)*sin(phi))^2 ) );
% B_CAN_Y = ( B_CAN_z  - (sin(theta)*sin(phi))/sqrt( 1-(sin(theta)*sin(phi))^2 ) * B_CAN_rho );
% 
% B_CAN_r = B_CAN_r + sin(theta)*sin(phi)*B_CAN_Y;
% B_CAN_theta = cos(theta)*sin(phi)*B_CAN_Y;
% B_CAN_phi = cos(phi)*B_CAN_Y;
%}
