function [FringingField, IsFirstIteration] = FringingField(k, r, ThetaMaxDeg, PhiMaxDeg, DeltaTheta, DeltaPhi, SystemParameters, IsFirstIteration)

    [NbPointsTotal, NbPointsInside,...
    N_theta, N_phi,...
    theta_vec, phi_vec] = GridDetails(ThetaMaxDeg, PhiMaxDeg, DeltaTheta*180/pi, DeltaPhi*180/pi);

    SurfaceGrid = reshape(r, [N_phi, N_theta]);
    
    FringingField = zeros(N_phi, N_theta);
    
    if IsFirstIteration
        
        FringingField = 0;
        
    else

    if k == 1
        n_Spherical = [ 1, -1/r(k)*(r(k+N_phi)-rss) / (2*delta_theta), -1/(r(k)*sin(theta_vec(k))) * (r(k+1)-r_eq(2)) / (2*delta_phi)];
        BDipole_Spherical = M*(1/r(k))^3 * [ -2*sin(theta_vec(k))*sin(phi_vec(k)), cos(theta_vec(k))*sin(phi_vec(k)), cos(phi_vec(k)) ];
        
        BfPrevious = FringingField(k, r, Max_deg_equ, Max_deg_phi_direction, delta_theta_deg, delta_phi_deg, SystemParameters, IsFirstIteration);        
        Btot = BDipole_Spherical + BfPrevious;
        CurrentDensity = - cross(n_Spherical, Btot);
        
        
        dBfNew = (1/(2*pi)) * CurrentDensity
        
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
        
        n_Spherical = [ 1, -1/r(k)*(r(k)pN-rss) / (2*delta_theta), -1/(r(k)*sin(theta_vec(k))) * (r_mer(2)-r(k)m1) / (2*delta_phi)];
%         n_Spherical = [ 1, -1/r(k)*(r(k)-rss) / (delta_theta), -1/(r(k)*sin(theta_vec(k))) * (r_mer(2)-r(k)m1) / (2*delta_phi)];
        BDipole_Spherical = M*(1/r(k))^3 * [ -2*sin(theta_vec(k))*sin(phi_vec(k)), cos(theta_vec(k))*sin(phi_vec(k)), cos(phi_vec(k)) ];

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

        n_Spherical = [ 1, -1/r(k)*(2*r(k)-r(k)mN-r(k)mN)/ (2*delta_theta), -1/(r(k)*sin(theta_vec(k))) * (r(k)p1-r_eq(N_theta+1)) / (2*delta_phi)];
%         n_Spherical = [ 1, -1/r(k)*(r(k)-r(k)mN)/ (2*delta_theta), -1/(r(k)*sin(theta_vec(k))) * (r(k)p1-r_eq(N_theta+1)) / (2*delta_phi)];
        BDipole_Spherical = M*(1/r(k))^3 * [ -2*sin(theta_vec(k))*sin(phi_vec(k)), cos(theta_vec(k))*sin(phi_vec(k)), cos(phi_vec(k)) ];

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
        
        n_Spherical = [ 1, -1/r(k)*(2*r(k)-r(k)mN-r(k)mN)/ (2*delta_theta), -1/(r(k)*sin(theta_vec(k))) * (r_mer(N_theta+1)-r(k)m1) / (2*delta_phi)];
%         n_Spherical = [ 1, -1/r(k)*(r(k)-r(k)mN)/ (delta_theta), -1/(r(k)*sin(theta_vec(k))) * (r_mer(N_theta+1)-r(k)m1) / (2*delta_phi)];
        BDipole_Spherical = M*(1/r(k))^3 * [ -2*sin(theta_vec(k))*sin(phi_vec(k)), cos(theta_vec(k))*sin(phi_vec(k)), cos(phi_vec(k)) ];

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

        n_Spherical = [ 1, -1/r(k)*(r(k)pN-rss) / (2*delta_theta), -1/(r(k)*sin(theta_vec(k))) * (r(k)p1-r(k)m1) / (2*delta_phi)];
%         n_Spherical = [ 1, -1/r(k)*(r(k)-rss) / (delta_theta), -1/(r(k)*sin(theta_vec(k))) * (r(k)p1-r(k)m1) / (2*delta_phi)];
        BDipole_Spherical = M*(1/r(k))^3 * [ -2*sin(theta_vec(k))*sin(phi_vec(k)), cos(theta_vec(k))*sin(phi_vec(k)), cos(phi_vec(k)) ];

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

        n_Spherical = [ 1, -1/r(k)*(r(k)pN-r(k)mN) / (2*delta_theta), -1/(r(k)*sin(theta_vec(k))) * (r(k)p1-r_eq((k-1)/N_phi+2)) / (2*delta_phi)];
%         n_Spherical = [ 1, -1/r(k)*(r(k)-r(k)mN) / (delta_theta), -1/(r(k)*sin(theta_vec(k))) * (r(k)p1-r_eq((k-1)/N_phi+2)) / (2*delta_phi)];
        BDipole_Spherical = M*(1/r(k))^3 * [ -2*sin(theta_vec(k))*sin(phi_vec(k)), cos(theta_vec(k))*sin(phi_vec(k)), cos(phi_vec(k)) ];

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
        
        n_Spherical = [ 1, -1/r(k)*(r(k)pN-r(k)mN) / (2*delta_theta), -1/(r(k)*sin(theta_vec(k))) * (r_mer(k/N_phi+1)-r(k)m1) / (2*delta_phi)];
%         n_Spherical = [ 1, -1/r(k)*(r(k)-r(k)mN) / (delta_theta), -1/(r(k)*sin(theta_vec(k))) * (r_mer(k/N_phi+1)-r(k)m1) / (2*delta_phi)];
    BDipole_Spherical = M*(1/r(k))^3 * [ -2*sin(theta_vec(k))*sin(phi_vec(k)), cos(theta_vec(k))*sin(phi_vec(k)), cos(phi_vec(k)) ];

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

        n_Spherical = [ 1, -1/r(k)*(2*r(k)-r(k)mN-r(k)mN) / (2*delta_theta), -1/(r(k)*sin(theta_vec(k))) * (r(k)p1-r(k)m1) / (2*delta_phi)];
%         n_Spherical = [ 1, -1/r(k)*(r(k)-r(k)mN) / (delta_theta), -1/(r(k)*sin(theta_vec(k))) * (r(k)p1-r(k)m1) / (2*delta_phi)];
        BDipole_Spherical = M*(1/r(k))^3 * [ -2*sin(theta_vec(k))*sin(phi_vec(k)), cos(theta_vec(k))*sin(phi_vec(k)), cos(phi_vec(k)) ];

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

        n_Spherical = [ 1, -1/r(k)*(r(k)pN-r(k)mN) / (2*delta_theta), -1/(r(k)*sin(theta_vec(k))) * (r(k)p1-r(k)m1) / (2*delta_phi)];
%         n_Spherical = [ 1, -1/r(k)*(r(k)-r(k)mN) / (delta_theta), -1/(r(k)*sin(theta_vec(k))) * (r(k)p1-r(k)m1) / (2*delta_phi)];
        BDipole_Spherical = M*(1/r(k))^3 * [ -2*sin(theta_vec(k))*sin(phi_vec(k)), cos(theta_vec(k))*sin(phi_vec(k)), cos(phi_vec(k)) ];

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
        
        
        
        
        
    end
    
    IsFirstIteration = "False";

end