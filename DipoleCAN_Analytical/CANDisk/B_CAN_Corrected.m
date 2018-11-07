function B_CAN_Corrected = B_CAN_Corrected(CAN_DiskParameters, r, theta, phi)
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

% Converting the CAN disk field components in MB64 Coordinate system
    if theta == 0 && phi == 0
        B_CAN_r = 0;
        B_CAN_theta = B_CAN_z;
        B_CAN_phi = 0;   
    elseif theta == pi/2 && phi == pi/2
        B_CAN_r = B_CAN_z;
        B_CAN_theta = 0;
        B_CAN_phi = 0;  
    elseif phi == 0
        B_CAN_r = 0;
        B_CAN_theta = 0;
        B_CAN_phi = B_CAN_z;   
    elseif phi == pi/2
        B_CAN_r = B_CAN_rho*cos(theta) + B_CAN_z*sin(theta);
        B_CAN_theta = -B_CAN_rho*sin(theta) + B_CAN_z*cos(theta);
        B_CAN_phi = 0;   
% 
%     elseif theta == pi/2 && phi == pi/2
%         B_CAN_r = B_CAN_z;
%         B_CAN_theta = 0;
%         B_CAN_phi = 0;   
% %     elseif theta <= 0*pi/180
% %         B_CAN_r = 0;
% %         B_CAN_theta = 0;
%         B_CAN_phi = B_CAN_z;
%         
%     elseif phi == 0
%         B_CAN_r = 0;
%         B_CAN_theta = B_CAN_z;
%         B_CAN_phi = 0;  
%         
%     elseif phi == pi/2
%         X_MB = 0;
%         Y_MB = B_CAN_z;
%         Z_MB = B_CAN_rho;
%         
% %         [B_CAN_phi, ~, B_CAN_r] = cart2sph(X_MB, Y_MB, Z_MB);
%         [B_CAN_phi, B_CAN_r] = cart2pol(Z_MB, Y_MB);
%         B_CAN_theta = 0;
% %         B_CAN_theta = atan2(sqrt(X_MB^2+ Y_MB^2), Z_MB);        
    else
        [X_CAN, Y_CAN, Z_CAN] = pol2cart(0, B_CAN_rho, B_CAN_z);
    
        X_MB = Y_CAN;
        Y_MB = Z_CAN;
        Z_MB = X_CAN;
    
        [B_CAN_phi, ~, B_CAN_r] = cart2sph(X_MB, Y_MB, Z_MB);

        B_CAN_theta = atan2(sqrt(X_MB^2+ Y_MB^2), Z_MB);
    end

    B_CAN_Corrected = [B_CAN_r; B_CAN_theta; B_CAN_phi];

return;    
