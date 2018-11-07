function [theta, r_eq] = Balance_Equator_CAN(theta_span, rSubSolarNose, beta)


%%%%% Along the equator, phi = 0, theta varies
options = odeset('RelTol',1e-6,'AbsTol',1e-6);
r_ini = rSubSolarNose;                                        % Initial "nose" distance, r0
dr_ini = 0;

[theta, r_eq] = ode15i( @(theta, r, dr) PressureBalanceCANDisk_drTheta(theta, 0, r, dr, beta), theta_span, r_ini, dr_ini, options );


end