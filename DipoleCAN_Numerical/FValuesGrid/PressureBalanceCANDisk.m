


function PressureBalanceCANDisk = PressureBalanceCANDisk (theta, phi, r, drdtheta, drdphi, beta, SystemParameters)

M = SystemParameters.M;


% Normal vector, Dipolar field and SW velocity vector
    if theta == 0
        n_vec = [ 1 ; 0 ; 0 ].';
    else
        n_vec = [ 1 ; (-1/r)*drdtheta ; -1/(r*(sin(theta)))*drdphi ].';
    end

    B = (M/r^3) .* [ -2*sin(theta)*sin(phi) ; cos(theta)*sin(phi) ; cos(phi) ].';

        P_sph2cart = [ sin(theta)*cos(phi)     sin(theta)*sin(phi)    cos(theta)    ;...
                       cos(theta)*cos(phi)     cos(theta)*sin(phi)    -sin(theta)   ;...   
                       -sin(phi)               cos(phi)               0             ...
                      ];

        n_Cartesian = n_vec * P_sph2cart;
        BDipole_Cartesian = B * P_sph2cart;
        v_Cartesian = [ 0, 0, -1 ]; 

    BDisk_Cartesian = B_CAN_Cartesian(SystemParameters, r, theta, phi).';
    Btot = BDipole_Cartesian + BDisk_Cartesian;
    
    PressureDifference = norm(cross(n_Cartesian, Btot))*sqrt(1+beta) + (1/2)*dot(n_Cartesian, v_Cartesian);
    MeanPressure = (1/2) * (norm(cross(n_Cartesian, Btot))*sqrt(1+beta) - (1/2)*dot(n_Cartesian,v_Cartesian));

%    PressureBalanceCANDisk = PressureDifference/MeanPressure;
    PressureBalanceCANDisk = PressureDifference;

% % System parameters
% M = SystemParameters.M;
%     
% % Normal vector, Dipolar field and SW velocity vector
%     if theta == 0
%         n_vec = [ 1 ; 0 ; 0 ];
%     else
%         n_vec = [ 1 ; (-1/r)*drdtheta ; -1/(r*(sin(theta)))*drdphi ];
%     end
%     
%     B = (M/r^3) .* [ -2*sin(theta)*sin(phi) ; cos(theta)*sin(phi) ; cos(phi) ];
%     v = [ -cos(theta) ; sin(theta) ; 0 ]; 
%     
% % CAN disk field in MB64 coordinate system
% if phi == 0 && theta == 0
%     B_Disk = B_CAN_Equator(SystemParameters, r, theta, phi);
% elseif phi == 0
%     B_Disk = B_CAN_Equator(SystemParameters, r, theta, phi);
% elseif phi == pi/2
%     B_Disk = B_CAN_Meridian(SystemParameters, r, theta, phi);
% end
% 
% % Total field in MB64 Coordinate system
%     Btot = B+B_Disk;
% 
% PressureBalanceCANDisk = norm(cross(n_vec, Btot))*sqrt(1+beta) + (1/2)*dot(n_vec,v);

end


 %{
%     
% % CAN disk field in CAN89 Coordinate system
%     z = r.*r0.*sin(theta).*sin(phi) / Rp ;
%     rho = r.*r0.*sqrt( 1-(sin(theta)*sin(phi))^2 ) / Rp ;
%     xdata = [z; rho];
%     B_CAN = diskfield_cyl(CAN_DiskParameters, xdata);
%     
%     B_CAN_z = B_CAN(1);
%     B_CAN_rho = B_CAN(2);
% 
% % CAN disk field in MB64 Coordinate system
%     if theta == pi/2 && phi == pi/2
%         B_CAN_r = B_CAN_z;
%         B_CAN_theta = 0;
%         B_CAN_phi = 0;   
%     elseif phi == 0
%         B_CAN_r = B_CAN_rho;
%         B_CAN_theta = 0;
%         B_CAN_phi = B_CAN_z;
%     else
%         [X_CAN, Y_CAN, Z_CAN] = pol2cart(0, B_CAN_rho, B_CAN_z);
%     
%         X_MB = Y_CAN;
%         Y_MB = Z_CAN;
%         Z_MB = X_CAN;
%     
%         [B_CAN_phi, ~, B_CAN_r] = cart2sph(X_MB, Y_MB, Z_MB);
% 
%         B_CAN_theta = atan2(sqrt(X_MB^2+ Y_MB^2), Z_MB);
%     end
% 
%     B_Disk = [B_CAN_r; B_CAN_theta; B_CAN_phi];

%}

