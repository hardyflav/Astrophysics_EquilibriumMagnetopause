function [theta, r_eq] = Balance_Equator(theta_span)

%%%%% System Data - Saturn
    mu0 = 4*pi*10^-7;                                 % Vaccum Permeability, N/A^2
    Bp = 20000*10^(-9);     % Equatorial field, T
    Rp = 60280*10^3;        % Planet radius, m
    
    Pswnpa = 0.02;
    
    Rsat = 60268e3;                                   % Saturn's equatorial radius, m
    rho_sw = ( 7*10^(6) ) * ( 1.6e-27 );              % Solar wind mass density, kg/m^3
    u_sw = 400e3;                                     % Solar wind velocity, in m/s

%%%%% Scaling coefficients
    b0 = sqrt(2*mu0*Pswnpa*10^(-9));    % Field scale
    r0 = (2*Bp*Rp^3/b0)^(1/3);          % Distance scale
    M = (Bp/b0)*(Rp/r0)^3;


%%%%% Along the equator, phi = 0, theta varies
options = odeset('RelTol',1e-10,'AbsTol',1e-10);
r_ini = 1;                                        % Initial "nose" distance, r0
dr_ini = 0;

[theta, r_eq] = ode15i( @(theta, r, dr) balance_eq_3D_theta(theta, 0, r, dr), theta_span, r_ini, dr_ini, options );


end