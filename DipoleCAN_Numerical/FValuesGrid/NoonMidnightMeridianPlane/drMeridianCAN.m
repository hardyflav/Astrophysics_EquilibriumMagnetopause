function drMeridianCAN = drMeridianCAN(theta, r, SystemParameters)

    CAN_DiskParameters = SystemParameters.CAN_DiskParameters;
    beta = SystemParameters.beta;
    M = SystemParameters.M;

    PhiMeridian = pi/2;
    B_Disk = B_CAN_Spherical(CAN_DiskParameters, r, theta, PhiMeridian);
        
    B_CAN_r = B_Disk(1);
    B_CAN_theta = B_Disk(2);

drMeridianCAN = (4.*B_CAN_r.^2.*(1+beta).*r.^6+16.*B_CAN_r.*(1+beta).*M.*r.^3.*sin( ...
  theta)+(16.*(1+beta).*M.^2+(-1).*r.^6).*sin(theta).^2).^(-1).*(r.* ...
  cos(theta).*(4.*B_CAN_r.*(1+beta).*M.*r.^3+(8.*(1+beta).*M.^2+r.^6) ...
  .*sin(theta))+(-2).*(2.*B_CAN_r.*B_CAN_theta.*(1+beta).*r.^7+4.* ...
  B_CAN_theta.*(1+beta).*M.*r.^4.*sin(theta)+((1+beta).*r.^8.*( ...
  B_CAN_theta.*r.^3.*sin(theta)+(-1).*cos(theta).*(B_CAN_r.*r.^3+3.*M.* ...
  sin(theta))).^2).^(1/2)));

end