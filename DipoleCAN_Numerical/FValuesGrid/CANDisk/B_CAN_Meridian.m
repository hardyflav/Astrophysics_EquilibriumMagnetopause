function B_CAN_Meridian = B_CAN_Meridian(SystemParameters, r, theta, phi)
%
% provides cylindrical components of field for axisymmetric ring current
% extending from a to b, with half-thickness D, and current parameter I
%
% input: x = [I a b D]
%        xdata = [z; rho] in 2 column format (phi is irrelevant)
%
% output: Bcyl_disk, contains both components sequentially
% ie
% Bz_disk = Bcyl_disk(1:length(Bcyl_disk)/2);
% Brho_disk= Bcyl_disk(length(Bcyl_disk)/2+1:length(Bcyl_disk));
%
% Giacomo Giampieri, October 2002


% %%% System parameters    
r0 = SystemParameters.r0;
Rp = SystemParameters.Rp;
CAN_DiskParameters = SystemParameters.CAN_DiskParameters;

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

% Converting the CAN disk field components in MB64 Coordinate system
    if theta == 0 && phi == 0
        B_CAN_r = 0;
        B_CAN_theta = B_CAN_z;
        B_CAN_phi = 0;   
    elseif theta == pi/2 && phi == pi/2
        B_CAN_r = B_CAN_z;
        B_CAN_theta = 0;
        B_CAN_phi = 0;   
    elseif phi == pi/2
        B_CAN_r = B_CAN_rho*cos(theta) + B_CAN_z*sin(theta);
        B_CAN_theta = -B_CAN_rho*sin(theta) + B_CAN_z*cos(theta);
        B_CAN_phi = 0;     
    end
    
    B_CAN_Meridian = [B_CAN_r; B_CAN_theta; B_CAN_phi];

end    
