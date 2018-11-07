function dBCANdr = dBCANdr(CAN_DiskParameters, r, theta, phi)

% System parameters
    Bp = 20000*10^(-9);     % Equatorial field, T
    Rp = 60280*10^3;        % Planet radius, m
    
    mu_0 = 4*pi*10^(-7);
    Pswnpa = 0.02;
    
    b0 = sqrt(2*mu_0*Pswnpa*10^(-9));   % Field scale
    r0 = (2*Bp*Rp^3/b0)^(1/3);          % Distance scale
    M = (Bp/b0)*(Rp/r0)^3;

% Converting the position from MB64 to CAN89, in Rp
    if theta == 0
        z = 0;
        rho = r *r0/Rp;
    elseif theta == pi/2 && phi == pi/2
        z = r *r0/Rp;
        rho = 0;
    else
         z = r.*sin(theta).*sin(phi)  *r0/Rp ;
        rho = r.*sqrt( 1-(sin(theta)*sin(phi))^2 )  *r0/Rp;
    end

% Computing the CAN disk field components in CAN89
    x = CAN_DiskParameters;
    xdata = [z; rho];
    B_CAN = diskfield_cyl(x, xdata);
    B_CAN_z = B_CAN(1);
    B_CAN_rho = B_CAN(2);

% Computing the derivative of the CAN disk field components in CAN89
dBCANdrho ?
dBCANdz ?

dBCANdr_CAN89 = (rho/r) * dBCANdrho + (z/r) * dBCANdz;





dBCANdr = [dBCANdr_r, dBCANdr_theta, dBCANdr_phi];

end
