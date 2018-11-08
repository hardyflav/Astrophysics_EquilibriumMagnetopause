function [PB_Cartesian_Numerical] = PB_Cartesian_Numerical( k, rk, rkp1, rkm1, rkpN, rkmN, rss, r_eq, r_mer, delta_theta_deg, delta_phi_deg, Max_deg_equ, Max_deg_phi_direction, SystemParameters )
    
    beta = SystemParameters.beta;
    M = SystemParameters.M;
    mu0I = SystemParameters.mu0I;
    
    [~, ~,...
    N_theta, N_phi,...
    theta_vec, phi_vec] = GridDetails(Max_deg_equ, Max_deg_phi_direction, delta_theta_deg, delta_phi_deg);

    delta_theta = delta_theta_deg*(pi/180);
    delta_phi = delta_phi_deg*(pi/180);
    
    if k == 1
        
        n_Spherical = [ 1, -1/rk*(rkpN-rss) / (2*delta_theta), -1/(rk*sin(theta_vec(k))) * (rkp1-r_eq(2)) / (2*delta_phi)];
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
        
        n_Spherical = [ 1, -1/rk*(rkpN-rss) / (2*delta_theta), -1/(rk*sin(theta_vec(k))) * (r_mer(2)-rkm1) / (2*delta_phi)];
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

        n_Spherical = [ 1, -1/rk*(rkpN-rss) / (2*delta_theta), -1/(rk*sin(theta_vec(k))) * (rkp1-rkm1) / (2*delta_phi)];
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
    
    if mu0I == 0
        Btot = BDipole_Cartesian;
    else
        BDisk_Cartesian = B_CAN_Cartesian(SystemParameters, rk, theta_vec(k), phi_vec(k)).';
        Btot = BDipole_Cartesian + BDisk_Cartesian;
    end
    
    PressureDifference = norm(cross(n_Cartesian, Btot))*sqrt(1+beta) + (1/2)*dot(n_Cartesian,v_Cartesian);
    MeanPressure = (1/2) * (norm(cross(n_Cartesian, Btot))*sqrt(1+beta) - (1/2)*dot(n_Cartesian,v_Cartesian));
%     PB_Cartesian_Numerical = (PressureDifference / MeanPressure);
    PB_Cartesian_Numerical = PressureDifference;
    
end
