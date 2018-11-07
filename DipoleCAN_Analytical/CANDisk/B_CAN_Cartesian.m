function B_CAN_Cartesian = B_CAN_Cartesian(SystemParameters, r, theta, phi)
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

% System parameters
    Rp = SystemParameters.Rp;
    r0 = SystemParameters.r0;
    CAN_DiskParameters = SystemParameters.CAN_DiskParameters;


% Converting the position from spherical MB64 to cartesian MB64, in Rp
    Azymuth = phi;
    Elevation = pi/2-theta;
    RadialPosition = r*r0/Rp;
    [X_MB, Y_MB, Z_MB] = sph2cart(Azymuth, Elevation, RadialPosition);
    
% Converting the position from cartesian MB64 to polar CAN89, in Rp
    [phi_CAN, rho_CAN, z_CAN] = cart2pol(Z_MB, X_MB, Y_MB);
    
% Computing the CAN disk field components in CAN89
    x = CAN_DiskParameters;
    xdata = [z_CAN; rho_CAN];
    B_CAN = diskfield_cyl(x, xdata);
    
    P_pol2cart = [ cos(phi_CAN)     sin(phi_CAN)    0    ;...
                  -sin(phi_CAN)     cos(phi_CAN)    0    ;...   
                  0                 0               1     ...
                 ];
    
    B_CAN_Cartesian = [B_CAN(2), 0, B_CAN(1)] *  P_pol2cart;
    B_CAN_X = B_CAN_Cartesian(1);
    B_CAN_Y = B_CAN_Cartesian(2);
    B_CAN_Z = B_CAN_Cartesian(3);    
    
% Converting the CAN disk field components in cartesian MB64
    BDisk_Z = B_CAN_X;
    BDisk_X = B_CAN_Y;
    BDisk_Y = B_CAN_Z;

    
    B_CAN_Cartesian = [BDisk_X; BDisk_Y; BDisk_Z];

return;    
