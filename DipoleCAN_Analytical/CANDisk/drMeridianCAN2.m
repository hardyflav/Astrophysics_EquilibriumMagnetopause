function drMeridianCAN2 = drMeridianCAN2(theta, r, beta, CAN_DiskParameters)

    Bp = 20000*10^(-9);     % Equatorial field, T
    Rp = 60280*10^3;        % Planet radius, m
    
    mu_0 = 4*pi*10^(-7);
    Pswnpa = 0.02;
    
    b0 = sqrt(2*mu_0*Pswnpa*10^(-9));   % Field scale
    r0 = (2*Bp*Rp^3/b0)^(1/3);          % Distance scale
    M = (Bp/b0)*(Rp/r0)^3;

    PhiMeridian = pi/2;
    B_Disk = B_CAN_Corrected(CAN_DiskParameters, r, theta, PhiMeridian);
    B_CAN_r = B_Disk(1);
    B_CAN_theta = B_Disk(2);


drMeridianCAN2 = (4.*B_CAN_r.^2.*(1+beta).*r.^6+16.*B_CAN_r.*(1+beta).*M.*r.^3.*sin( ...
  theta)+(16.*(1+beta).*M.^2+(-1).*r.^6).*sin(theta).^2).^(-1).*(( ...
  -4).*B_CAN_r.*B_CAN_theta.*(1+beta).*r.^7+(-8).*B_CAN_theta.*(1+beta).* ...
  M.*r.^4.*sin(theta)+r.*cos(theta).*(4.*B_CAN_r.*(1+beta).*M.*r.^3+( ...
  8.*(1+beta).*M.^2+r.^6).*sin(theta))+2.*((1+beta).*r.^8.*( ...
  B_CAN_theta.*r.^3.*sin(theta)+(-1).*cos(theta).*(B_CAN_r.*r.^3+3.*M.* ...
  sin(theta))).^2).^(1/2));

end