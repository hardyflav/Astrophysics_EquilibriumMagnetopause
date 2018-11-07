function [PB_Cartesian] = PB_Cartesian( k, r, rss, r_eq, r_mer, delta_theta_deg, delta_phi_deg, Max_deg_equ, Max_deg_phi_direction, beta, SystemParameters )
    

%     M = SystemParameters.M;
%     
    [~, ~,...
    N_theta, N_phi,...
    ~, ~] = GridDetails(Max_deg_equ, Max_deg_phi_direction, delta_theta_deg, delta_phi_deg);

%     delta_theta = delta_theta_deg*(pi/180);
%     delta_phi = delta_phi_deg*(pi/180);
    
    if k == 1
        
        PB_Cartesian = PB_Cartesian_Numerical( k, r(k), r(k+1), r_eq(2), r(k+N_phi), rss, rss, r_eq, r_mer, delta_theta_deg, delta_phi_deg, Max_deg_equ, Max_deg_phi_direction, SystemParameters );
% 
%         n_Spherical = [ 1, -1/r(k)*(r(k+N_phi)-rss) / (2*delta_theta), -1/(r(k)*sin(theta_vec(k))) * (r(k+1)-r_eq(2)) / (2*delta_phi)];
% %         n_Spherical = [ 1, -1/r(k)*(r(k)-rss) / (delta_theta), -1/(r(k)*sin(theta_vec(k))) * (r(k)-r_eq(2)) / (delta_phi)];
%         BDipole_Spherical = M*(1/r(k))^3 * [ -2*sin(theta_vec(k))*sin(phi_vec(k)), cos(theta_vec(k))*sin(phi_vec(k)), cos(phi_vec(k)) ];
%         
%         theta = theta_vec(k);
%         phi = phi_vec(k);
%         
%         P_sph2cart = [ sin(theta)*cos(phi)     sin(theta)*sin(phi)    cos(theta)    ;...
%                        cos(theta)*cos(phi)     cos(theta)*sin(phi)    -sin(theta)   ;...   
%                        -sin(phi)               cos(phi)               0             ...
%                       ];
%         
%         n_Cartesian = n_Spherical * P_sph2cart;
%         BDipole_Cartesian = BDipole_Spherical * P_sph2cart;
%         v_Cartesian = [ 0, 0, -1 ];        

        
    elseif k == N_phi
        
        PB_Cartesian = PB_Cartesian_Numerical( k, r(k), r_mer(2), r(k-1), r(k+N_phi), rss, rss, r_eq, r_mer, delta_theta_deg, delta_phi_deg, Max_deg_equ, Max_deg_phi_direction, SystemParameters );
%         
%         n_Spherical = [ 1, -1/r(k)*(r(k+N_phi)-rss) / (2*delta_theta), -1/(r(k)*sin(theta_vec(k))) * (r_mer(2)-r(k-1)) / (2*delta_phi)];
% %         n_Spherical = [ 1, -1/r(k)*(r(k)-rss) / (delta_theta), -1/(r(k)*sin(theta_vec(k))) * (r(k)-r(k-1)) / (delta_phi)];
%         BDipole_Spherical = M*(1/r(k))^3 * [ -2*sin(theta_vec(k))*sin(phi_vec(k)), cos(theta_vec(k))*sin(phi_vec(k)), cos(phi_vec(k)) ];
% 
%         theta = theta_vec(k);
%         phi = phi_vec(k);
%         
%         P_sph2cart = [ sin(theta)*cos(phi)     sin(theta)*sin(phi)    cos(theta)    ;...
%                        cos(theta)*cos(phi)     cos(theta)*sin(phi)    -sin(theta)   ;...   
%                        -sin(phi)               cos(phi)               0             ...
%                       ];
%         
%         n_Cartesian = n_Spherical * P_sph2cart;
%         BDipole_Cartesian = BDipole_Spherical * P_sph2cart;
%         v_Cartesian = [ 0, 0, -1 ];        

    elseif k == N_phi*N_theta-N_phi+1

        PB_Cartesian = PB_Cartesian_Numerical( k, r(k), r(k+1), r_eq(N_theta+1), 2*r(k)-r(k-N_phi), r(k-N_phi), rss, r_eq, r_mer, delta_theta_deg, delta_phi_deg, Max_deg_equ, Max_deg_phi_direction, SystemParameters );
%         
%         n_Spherical = [ 1, -1/r(k)*(2*r(k)-r(k-N_phi)-r(k-N_phi))/ (2*delta_theta), -1/(r(k)*sin(theta_vec(k))) * (r(k+1)-r_eq(N_theta+1)) / (2*delta_phi)];
% %         n_Spherical = [ 1, -1/r(k)*(r(k)-r(k-N_phi))/ (delta_theta), -1/(r(k)*sin(theta_vec(k))) * (r(k)-r_eq(N_theta+1)) / (delta_phi)];
%         BDipole_Spherical = M*(1/r(k))^3 * [ -2*sin(theta_vec(k))*sin(phi_vec(k)), cos(theta_vec(k))*sin(phi_vec(k)), cos(phi_vec(k)) ];
% 
%         theta = theta_vec(k);
%         phi = phi_vec(k);
%         
%         P_sph2cart = [ sin(theta)*cos(phi)     sin(theta)*sin(phi)    cos(theta)    ;...
%                        cos(theta)*cos(phi)     cos(theta)*sin(phi)    -sin(theta)   ;...   
%                        -sin(phi)               cos(phi)               0             ...
%                       ];
%         
%         n_Cartesian = n_Spherical * P_sph2cart;
%         BDipole_Cartesian = BDipole_Spherical * P_sph2cart;
%         v_Cartesian = [ 0, 0, -1 ];        
    
    elseif k == N_phi*N_theta
        
        PB_Cartesian = PB_Cartesian_Numerical( k, r(k), r_mer(N_theta+1), r(k-1), 2*r(k)-r(k-N_phi), r(k-N_phi), rss, r_eq, r_mer, delta_theta_deg, delta_phi_deg, Max_deg_equ, Max_deg_phi_direction, SystemParameters );
        
%         n_Spherical = [ 1, -1/r(k)*(2*r(k)-r(k-N_phi)-r(k-N_phi))/ (2*delta_theta), -1/(r(k)*sin(theta_vec(k))) * (r_mer(N_theta+1)-r(k-1)) / (2*delta_phi)];
% %         n_Spherical = [ 1, -1/r(k)*(2*r(k)-r(k)-r(k-N_phi))/ (delta_theta), -1/(r(k)*sin(theta_vec(k))) * (r(k)-r(k-1)) / (delta_phi)];
%         BDipole_Spherical = M*(1/r(k))^3 * [ -2*sin(theta_vec(k))*sin(phi_vec(k)), cos(theta_vec(k))*sin(phi_vec(k)), cos(phi_vec(k)) ];
% 
%         theta = theta_vec(k);
%         phi = phi_vec(k);
%         
%         P_sph2cart = [ sin(theta)*cos(phi)     sin(theta)*sin(phi)    cos(theta)    ;...
%                        cos(theta)*cos(phi)     cos(theta)*sin(phi)    -sin(theta)   ;...   
%                        -sin(phi)               cos(phi)               0             ...
%                       ];
%         
%         n_Cartesian = n_Spherical * P_sph2cart;
%         BDipole_Cartesian = BDipole_Spherical * P_sph2cart;
%         v_Cartesian = [ 0, 0, -1 ];   
         
    elseif k>1 && k<N_phi   %%% left boundary

        PB_Cartesian = PB_Cartesian_Numerical( k, r(k), r(k+1), r(k-1), r(k+N_phi), rss, rss, r_eq, r_mer, delta_theta_deg, delta_phi_deg, Max_deg_equ, Max_deg_phi_direction, SystemParameters );
        
%         n_Spherical = [ 1, -1/r(k)*(r(k+N_phi)-rss) / (2*delta_theta), -1/(r(k)*sin(theta_vec(k))) * (r(k+1)-r(k-1)) / (2*delta_phi)];
% %         n_Spherical = [ 1, -1/r(k)*(r(k)-rss) / (delta_theta), -1/(r(k)*sin(theta_vec(k))) * (r(k)-r(k-1)) / (delta_phi)];
%         BDipole_Spherical = M*(1/r(k))^3 * [ -2*sin(theta_vec(k))*sin(phi_vec(k)), cos(theta_vec(k))*sin(phi_vec(k)), cos(phi_vec(k)) ];
% 
%         theta = theta_vec(k);
%         phi = phi_vec(k);
%         
%         P_sph2cart = [ sin(theta)*cos(phi)     sin(theta)*sin(phi)    cos(theta)    ;...
%                        cos(theta)*cos(phi)     cos(theta)*sin(phi)    -sin(theta)   ;...   
%                        -sin(phi)               cos(phi)               0             ...
%                       ];
%         
%         n_Cartesian = n_Spherical * P_sph2cart;
%         BDipole_Cartesian = BDipole_Spherical * P_sph2cart;
%         v_Cartesian = [ 0, 0, -1 ];  
    
    elseif k>1 && mod(k,N_phi)==1 && k<N_phi*N_theta-N_phi+1    %%% bottom boundary
        
        PB_Cartesian = PB_Cartesian_Numerical( k, r(k), r(k+1), r_eq((k-1)/N_phi+2), r(k+N_phi), r(k-N_phi), rss, r_eq, r_mer, delta_theta_deg, delta_phi_deg, Max_deg_equ, Max_deg_phi_direction, SystemParameters );

%         n_Spherical = [ 1, -1/r(k)*(r(k+N_phi)-r(k-N_phi)) / (2*delta_theta), -1/(r(k)*sin(theta_vec(k))) * (r(k+1)-r_eq((k-1)/N_phi+2)) / (2*delta_phi)];
% %         n_Spherical = [ 1, -1/r(k)*(r(k)-r(k-N_phi)) / (delta_theta), -1/(r(k)*sin(theta_vec(k))) * (r(k)-r_eq((k-1)/N_phi+2)) / (delta_phi)];
%         BDipole_Spherical = M*(1/r(k))^3 * [ -2*sin(theta_vec(k))*sin(phi_vec(k)), cos(theta_vec(k))*sin(phi_vec(k)), cos(phi_vec(k)) ];
% 
%         theta = theta_vec(k);
%         phi = phi_vec(k);
%         
%         P_sph2cart = [ sin(theta)*cos(phi)     sin(theta)*sin(phi)    cos(theta)    ;...
%                        cos(theta)*cos(phi)     cos(theta)*sin(phi)    -sin(theta)   ;...   
%                        -sin(phi)               cos(phi)               0             ...
%                       ];
%         
%         n_Cartesian = n_Spherical * P_sph2cart;
%         BDipole_Cartesian = BDipole_Spherical * P_sph2cart;
%         v_Cartesian = [ 0, 0, -1 ];  
        
    elseif k>N_phi && mod(k,N_phi)==0 && k<N_phi*N_theta  %%% upper boundary
        
        PB_Cartesian = PB_Cartesian_Numerical( k, r(k), r_mer(k/N_phi+1), r(k-1), r(k+N_phi), r(k-N_phi), rss, r_eq, r_mer, delta_theta_deg, delta_phi_deg, Max_deg_equ, Max_deg_phi_direction, SystemParameters );
        
%         n_Spherical = [ 1, -1/r(k)*(r(k+N_phi)-r(k-N_phi)) / (2*delta_theta), -1/(r(k)*sin(theta_vec(k))) * (r_mer(k/N_phi+1)-r(k-1)) / (2*delta_phi)];
% %         n_Spherical = [ 1, -1/r(k)*(r(k)-r(k-N_phi)) / (delta_theta), -1/(r(k)*sin(theta_vec(k))) * (r(k)-r(k-1)) / (delta_phi)];
%         BDipole_Spherical = M*(1/r(k))^3 * [ -2*sin(theta_vec(k))*sin(phi_vec(k)), cos(theta_vec(k))*sin(phi_vec(k)), cos(phi_vec(k)) ];
% 
%         theta = theta_vec(k);
%         phi = phi_vec(k);
%         
%         P_sph2cart = [ sin(theta)*cos(phi)     sin(theta)*sin(phi)    cos(theta)    ;...
%                        cos(theta)*cos(phi)     cos(theta)*sin(phi)    -sin(theta)   ;...   
%                        -sin(phi)               cos(phi)               0             ...
%                       ];
%         
%         n_Cartesian = n_Spherical * P_sph2cart;
%         BDipole_Cartesian = BDipole_Spherical * P_sph2cart;
%         v_Cartesian = [ 0, 0, -1 ];  
    
    elseif k>N_phi*N_theta-N_phi+1 && k<N_phi*N_theta  %%% right boundary
        
        PB_Cartesian = PB_Cartesian_Numerical( k, r(k), r(k+1), r(k-1), 2*r(k)-r(k-N_phi), r(k-N_phi), rss, r_eq, r_mer, delta_theta_deg, delta_phi_deg, Max_deg_equ, Max_deg_phi_direction, SystemParameters );

%         n_Spherical = [ 1, -1/r(k)*(2*r(k)-r(k-N_phi)-r(k-N_phi)) / (2*delta_theta), -1/(r(k)*sin(theta_vec(k))) * (r(k+1)-r(k-1)) / (2*delta_phi)];
% %         n_Spherical = [ 1, -1/r(k)*(r(k)-r(k-N_phi)) / (delta_theta), -1/(r(k)*sin(theta_vec(k))) * (r(k)-r(k-1)) / (delta_phi)];
%         BDipole_Spherical = M*(1/r(k))^3 * [ -2*sin(theta_vec(k))*sin(phi_vec(k)), cos(theta_vec(k))*sin(phi_vec(k)), cos(phi_vec(k)) ];
% 
%         theta = theta_vec(k);
%         phi = phi_vec(k);
%         
%         P_sph2cart = [ sin(theta)*cos(phi)     sin(theta)*sin(phi)    cos(theta)    ;...
%                        cos(theta)*cos(phi)     cos(theta)*sin(phi)    -sin(theta)   ;...   
%                        -sin(phi)               cos(phi)               0             ...
%                       ];
%         
%         n_Cartesian = n_Spherical * P_sph2cart;
%         BDipole_Cartesian = BDipole_Spherical * P_sph2cart;
%         v_Cartesian = [ 0, 0, -1 ];  
    
    else

        PB_Cartesian = PB_Cartesian_Numerical( k, r(k), r(k+1), r(k-1), r(k+N_phi), r(k-N_phi), rss, r_eq, r_mer, delta_theta_deg, delta_phi_deg, Max_deg_equ, Max_deg_phi_direction, SystemParameters );

%         n_Spherical = [ 1, -1/r(k)*(r(k+N_phi)-r(k-N_phi)) / (2*delta_theta), -1/(r(k)*sin(theta_vec(k))) * (r(k+1)-r(k-1)) / (2*delta_phi)];
% %         n_Spherical = [ 1, -1/r(k)*(r(k)-r(k-N_phi)) / (delta_theta), -1/(r(k)*sin(theta_vec(k))) * (r(k)-r(k-1)) / (delta_phi)];
%         BDipole_Spherical = M*(1/r(k))^3 * [ -2*sin(theta_vec(k))*sin(phi_vec(k)), cos(theta_vec(k))*sin(phi_vec(k)), cos(phi_vec(k)) ];
% 
%         theta = theta_vec(k);
%         phi = phi_vec(k);
%         
%         P_sph2cart = [ sin(theta)*cos(phi)     sin(theta)*sin(phi)    cos(theta)    ;...
%                        cos(theta)*cos(phi)     cos(theta)*sin(phi)    -sin(theta)   ;...   
%                        -sin(phi)               cos(phi)               0             ...
%                       ];
%         
%         n_Cartesian = n_Spherical * P_sph2cart;
%         BDipole_Cartesian = BDipole_Spherical * P_sph2cart;
%         v_Cartesian = [ 0, 0, -1 ];  
        
    end
    
%     BDisk_Cartesian = B_CAN_Cartesian(SystemParameters, r(k), theta_vec(k), phi_vec(k)).';
%     Btot = BDipole_Cartesian + BDisk_Cartesian;
%     
%     PressureDifference = norm(cross(n_Cartesian, Btot))*sqrt(1+beta) + (1/2)*dot(n_Cartesian,v_Cartesian);
%     MeanPressure = (1/2) * (norm(cross(n_Cartesian, Btot))*sqrt(1+beta) - (1/2)*dot(n_Cartesian,v_Cartesian));
% %     PB_Cartesian = PressureDifference / MeanPressure;
%     PB_Cartesian = PressureDifference;

end
