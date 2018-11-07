function [PB_Cartesian_Numerical_Bottom] = PB_Cartesian_Numerical_Bottom( k, rk, rkp1, rkm1, rkpN, rkmN, rLeft, r_eq, r_mer, delta_theta_deg, delta_phi_deg, Max_deg_equ, Max_deg_phi_direction, ThetaLeftBoundaryDeg, SystemParameters )
    
    beta = SystemParameters.beta;
    M = SystemParameters.M;
    
    Nb_points_theta = size(r_eq, 1);    % Number of points on each direction of the grid, from ThetaLeftBoundaryDeg to Max_deg_equ
    N_theta = Nb_points_theta-2; % Number of points on theta-direction of the grid, from 1 to Max_deg_equ-1

    Nb_points_phi= size(rLeft, 1);    % Number of points on each direction of the grid, from 0 to Max_deg_equ
    N_phi = Nb_points_phi-2; % Number of points on theta-direction of the grid, from 1 to Max_deg_equ-1

    delta_theta = delta_theta_deg*(pi/180);
    delta_phi = delta_phi_deg*(pi/180);
    
    ThetaMinDeg = ThetaLeftBoundaryDeg + delta_theta_deg;
    Max_deg_theta = Max_deg_equ-delta_theta_deg;
    Max_deg_phi= Max_deg_phi_direction-delta_phi_deg;
        
    theta_vec_1 = repmat(ThetaMinDeg:delta_theta_deg:Max_deg_theta,[N_phi 1]);
    theta_vec = theta_vec_1(:)*pi/180;

    phi_vec_1 = repmat(delta_phi_deg:delta_phi_deg:Max_deg_phi,[1 N_theta]).';
    phi_vec = phi_vec_1(:)*pi/180;
    
    rLeftBoundary = rLeft(2:end-1);

    if k == 1
        
        n_Spherical = [ 1, -1/rk*(rkpN-rLeftBoundary(1)) / (2*delta_theta), -1/(rk*sin(theta_vec(k))) * (rkp1-r_eq(2)) / (2*delta_phi)];
%         n_Spherical = [ 1, -1/rk*(rk-rLeftBoundary(1)) / (delta_theta), -1/(rk*sin(theta_vec(k))) * (rk-r_eq(2)) / (delta_phi)];
%         n_Spherical = [ 1, -1/rk*(rk-rss) / (delta_theta), -1/(rk*sin(theta_vec(k))) * (rkp1-r_eq(2)) / (2*delta_phi)];
        BDipole_Spherical = M*(1/rk)^3 * [ -2*sin(theta_vec(k))*sin(phi_vec(k)), cos(theta_vec(k))*sin(phi_vec(k)), cos(phi_vec(k)) ];
        
        theta = theta_vec(k);
        phi = phi_vec(k);
        
        P_sph2cart = [ sin(theta)*cos(phi)     sin(theta)*sin(phi)    cos(theta)    ;...
                       cos(theta)*cos(phi)     cos(theta)*sin(phi)    -sin(theta)   ;...   
                       -sin(phi)               cos(phi)               0             ...
                      ];
        
        n_Cartesian = n_Spherical * P_sph2cart;
        BDipole_Cartesian = BDipole_Spherical * P_sph2cart;
        v_Cartesian = [ 0, 0, -1 ];        

        
    elseif k == N_phi
        
        n_Spherical = [ 1, -1/rk*(rkpN-rLeftBoundary(end)) / (2*delta_theta), -1/(rk*sin(theta_vec(k))) * (r_mer(2)-rkm1) / (2*delta_phi)];
%         n_Spherical = [ 1, -1/rk*(rk-rLeftBoundary(end)) / (delta_theta), -1/(rk*sin(theta_vec(k))) * (rk-rkm1) / (delta_phi)];
%         n_Spherical = [ 1, -1/rk*(rk-rss) / (delta_theta), -1/(rk*sin(theta_vec(k))) * (r_mer(2)-rkm1) / (2*delta_phi)];
        BDipole_Spherical = M*(1/rk)^3 * [ -2*sin(theta_vec(k))*sin(phi_vec(k)), cos(theta_vec(k))*sin(phi_vec(k)), cos(phi_vec(k)) ];

        theta = theta_vec(k);
        phi = phi_vec(k);
        
        P_sph2cart = [ sin(theta)*cos(phi)     sin(theta)*sin(phi)    cos(theta)    ;...
                       cos(theta)*cos(phi)     cos(theta)*sin(phi)    -sin(theta)   ;...   
                       -sin(phi)               cos(phi)               0             ...
                      ];
        
        n_Cartesian = n_Spherical * P_sph2cart;
        BDipole_Cartesian = BDipole_Spherical * P_sph2cart;
        v_Cartesian = [ 0, 0, -1 ];        

    elseif k == N_phi*N_theta-N_phi+1

        n_Spherical = [ 1, -1/rk*(2*rk-rkmN-rkmN)/ (2*delta_theta), -1/(rk*sin(theta_vec(k))) * (rkp1-r_eq(N_theta+1)) / (2*delta_phi)];
%         n_Spherical = [ 1, -1/rk*(rk-rkmN)/ (delta_theta), -1/(rk*sin(theta_vec(k))) * (rk-r_eq(N_theta+1)) / (delta_phi)];
%         n_Spherical = [ 1, -1/rk*(rk-rkmN)/ (2*delta_theta), -1/(rk*sin(theta_vec(k))) * (rkp1-r_eq(N_theta+1)) / (2*delta_phi)];
        BDipole_Spherical = M*(1/rk)^3 * [ -2*sin(theta_vec(k))*sin(phi_vec(k)), cos(theta_vec(k))*sin(phi_vec(k)), cos(phi_vec(k)) ];

        theta = theta_vec(k);
        phi = phi_vec(k);
        
        P_sph2cart = [ sin(theta)*cos(phi)     sin(theta)*sin(phi)    cos(theta)    ;...
                       cos(theta)*cos(phi)     cos(theta)*sin(phi)    -sin(theta)   ;...   
                       -sin(phi)               cos(phi)               0             ...
                      ];
        
        n_Cartesian = n_Spherical * P_sph2cart;
        BDipole_Cartesian = BDipole_Spherical * P_sph2cart;
        v_Cartesian = [ 0, 0, -1 ];        
    
    elseif k == N_phi*N_theta
        
        n_Spherical = [ 1, -1/rk*(2*rk-rkmN-rkmN)/ (2*delta_theta), -1/(rk*sin(theta_vec(k))) * (r_mer(N_theta+1)-rkm1) / (2*delta_phi)];
%         n_Spherical = [ 1, -1/rk*(rk-rkmN)/ (delta_theta), -1/(rk*sin(theta_vec(k))) * (rk-rkm1) / (delta_phi)];
        %         n_Spherical = [ 1, -1/rk*(rk-rkmN)/ (delta_theta), -1/(rk*sin(theta_vec(k))) * (r_mer(N_theta+1)-rkm1) / (2*delta_phi)];
        BDipole_Spherical = M*(1/rk)^3 * [ -2*sin(theta_vec(k))*sin(phi_vec(k)), cos(theta_vec(k))*sin(phi_vec(k)), cos(phi_vec(k)) ];

        theta = theta_vec(k);
        phi = phi_vec(k);
        
        P_sph2cart = [ sin(theta)*cos(phi)     sin(theta)*sin(phi)    cos(theta)    ;...
                       cos(theta)*cos(phi)     cos(theta)*sin(phi)    -sin(theta)   ;...   
                       -sin(phi)               cos(phi)               0             ...
                      ];
        
        n_Cartesian = n_Spherical * P_sph2cart;
        BDipole_Cartesian = BDipole_Spherical * P_sph2cart;
        v_Cartesian = [ 0, 0, -1 ];   
        
    elseif k>1 && k<N_phi   %%% left boundary

        n_Spherical = [ 1, -1/rk*(rkpN-rLeftBoundary(k)) / (2*delta_theta), -1/(rk*sin(theta_vec(k))) * (rkp1-rkm1) / (2*delta_phi)];
%         n_Spherical = [ 1, -1/rk*(rk-rLeftBoundary(k)) / (delta_theta), -1/(rk*sin(theta_vec(k))) * (rk-rkm1) / (delta_phi)];
%         n_Spherical = [ 1, -1/rk*(rk-rss) / (delta_theta), -1/(rk*sin(theta_vec(k))) * (rkp1-rkm1) / (2*delta_phi)];
        BDipole_Spherical = M*(1/rk)^3 * [ -2*sin(theta_vec(k))*sin(phi_vec(k)), cos(theta_vec(k))*sin(phi_vec(k)), cos(phi_vec(k)) ];

        theta = theta_vec(k);
        phi = phi_vec(k);
        
        P_sph2cart = [ sin(theta)*cos(phi)     sin(theta)*sin(phi)    cos(theta)    ;...
                       cos(theta)*cos(phi)     cos(theta)*sin(phi)    -sin(theta)   ;...   
                       -sin(phi)               cos(phi)               0             ...
                      ];
        
        n_Cartesian = n_Spherical * P_sph2cart;
        BDipole_Cartesian = BDipole_Spherical * P_sph2cart;
        v_Cartesian = [ 0, 0, -1 ];  
    
    elseif k>1 && mod(k,N_phi)==1 && k<N_phi*N_theta-N_phi+1    %%% bottom boundary

        n_Spherical = [ 1, -1/rk*(rkpN-rkmN) / (2*delta_theta), -1/(rk*sin(theta_vec(k))) * (rkp1-r_eq((k-1)/N_phi+2)) / (2*delta_phi)];
%         n_Spherical = [ 1, -1/rk*(rk-rkmN) / (delta_theta), -1/(rk*sin(theta_vec(k))) * (rk-r_eq((k-1)/N_phi+2)) / (delta_phi)];
%         n_Spherical = [ 1, -1/rk*(rk-rkmN) / (delta_theta), -1/(rk*sin(theta_vec(k))) * (rkp1-r_eq((k-1)/N_phi+2)) / (2*delta_phi)];
        BDipole_Spherical = M*(1/rk)^3 * [ -2*sin(theta_vec(k))*sin(phi_vec(k)), cos(theta_vec(k))*sin(phi_vec(k)), cos(phi_vec(k)) ];

        theta = theta_vec(k);
        phi = phi_vec(k);
        
        P_sph2cart = [ sin(theta)*cos(phi)     sin(theta)*sin(phi)    cos(theta)    ;...
                       cos(theta)*cos(phi)     cos(theta)*sin(phi)    -sin(theta)   ;...   
                       -sin(phi)               cos(phi)               0             ...
                      ];
        
        n_Cartesian = n_Spherical * P_sph2cart;
        BDipole_Cartesian = BDipole_Spherical * P_sph2cart;
        v_Cartesian = [ 0, 0, -1 ];  
        
    elseif k>N_phi && mod(k,N_phi)==0 && k<N_phi*N_theta  %%% upper boundary
        
        n_Spherical = [ 1, -1/rk*(rkpN-rkmN) / (2*delta_theta), -1/(rk*sin(theta_vec(k))) * (r_mer(k/N_phi+1)-rkm1) / (2*delta_phi)];
%         n_Spherical = [ 1, -1/rk*(rk-rkmN) / (delta_theta), -1/(rk*sin(theta_vec(k))) * (rk-rkm1) / (delta_phi)];
%         n_Spherical = [ 1, -1/rk*(rk-rkmN) / (delta_theta), -1/(rk*sin(theta_vec(k))) * (r_mer(k/N_phi+1)-rkm1) / (2*delta_phi)];
    BDipole_Spherical = M*(1/rk)^3 * [ -2*sin(theta_vec(k))*sin(phi_vec(k)), cos(theta_vec(k))*sin(phi_vec(k)), cos(phi_vec(k)) ];

        theta = theta_vec(k);
        phi = phi_vec(k);
        
        P_sph2cart = [ sin(theta)*cos(phi)     sin(theta)*sin(phi)    cos(theta)    ;...
                       cos(theta)*cos(phi)     cos(theta)*sin(phi)    -sin(theta)   ;...   
                       -sin(phi)               cos(phi)               0             ...
                      ];
        
        n_Cartesian = n_Spherical * P_sph2cart;
        BDipole_Cartesian = BDipole_Spherical * P_sph2cart;
        v_Cartesian = [ 0, 0, -1 ];  
    
    elseif k>N_phi*N_theta-N_phi+1 && k<N_phi*N_theta  %%% right boundary

        n_Spherical = [ 1, -1/rk*(2*rk-rkmN-rkmN) / (2*delta_theta), -1/(rk*sin(theta_vec(k))) * (rkp1-rkm1) / (2*delta_phi)];
%         n_Spherical = [ 1, -1/rk*(rk-rkmN) / (delta_theta), -1/(rk*sin(theta_vec(k))) * (rk-rkm1) / (delta_phi)];
%         n_Spherical = [ 1, -1/rk*(rk-rkmN) / (delta_theta), -1/(rk*sin(theta_vec(k))) * (rkp1-rkm1) / (2*delta_phi)];
        BDipole_Spherical = M*(1/rk)^3 * [ -2*sin(theta_vec(k))*sin(phi_vec(k)), cos(theta_vec(k))*sin(phi_vec(k)), cos(phi_vec(k)) ];

        theta = theta_vec(k);
        phi = phi_vec(k);
        
        P_sph2cart = [ sin(theta)*cos(phi)     sin(theta)*sin(phi)    cos(theta)    ;...
                       cos(theta)*cos(phi)     cos(theta)*sin(phi)    -sin(theta)   ;...   
                       -sin(phi)               cos(phi)               0             ...
                      ];
        
        n_Cartesian = n_Spherical * P_sph2cart;
        BDipole_Cartesian = BDipole_Spherical * P_sph2cart;
        v_Cartesian = [ 0, 0, -1 ];  
    
    else

        n_Spherical = [ 1, -1/rk*(rkpN-rkmN) / (2*delta_theta), -1/(rk*sin(theta_vec(k))) * (rkp1-rkm1) / (2*delta_phi)];
%         n_Spherical = [ 1, -1/rk*(rk-rkmN) / (delta_theta), -1/(rk*sin(theta_vec(k))) * (rk-rkm1) / (delta_phi)];
%         n_Spherical = [ 1, -1/rk*(rk-rkmN) / (delta_theta), -1/(rk*sin(theta_vec(k))) * (rkp1-rkm1) / (2*delta_phi)];
        BDipole_Spherical = M*(1/rk)^3 * [ -2*sin(theta_vec(k))*sin(phi_vec(k)), cos(theta_vec(k))*sin(phi_vec(k)), cos(phi_vec(k)) ];

        theta = theta_vec(k);
        phi = phi_vec(k);
        
        P_sph2cart = [ sin(theta)*cos(phi)     sin(theta)*sin(phi)    cos(theta)    ;...
                       cos(theta)*cos(phi)     cos(theta)*sin(phi)    -sin(theta)   ;...   
                       -sin(phi)               cos(phi)               0             ...
                      ];
        
        n_Cartesian = n_Spherical * P_sph2cart;
        BDipole_Cartesian = BDipole_Spherical * P_sph2cart;
        v_Cartesian = [ 0, 0, -1 ];  
        
    end
    
    BDisk_Cartesian = B_CAN_Cartesian(SystemParameters, rk, theta_vec(k), phi_vec(k)).';
    Btot = BDipole_Cartesian + BDisk_Cartesian;
    
    PressureDifference = norm(cross(n_Cartesian, Btot))*sqrt(1+beta) + (1/2)*dot(n_Cartesian,v_Cartesian);
    MeanPressure = (1/2) * (norm(cross(n_Cartesian, Btot))*sqrt(1+beta) - (1/2)*dot(n_Cartesian,v_Cartesian));
%     PB_Cartesian_Numerical_Bottom = (PressureDifference / MeanPressure);
    PB_Cartesian_Numerical_Bottom = PressureDifference;
    
end
