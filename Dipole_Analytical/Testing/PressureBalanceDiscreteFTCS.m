function [PressureBalanceDiscreteFTCS] = PressureBalanceDiscreteFTCS( k, r, rss, r_eq, r_mer, delta_theta_deg, delta_phi_deg, Max_deg_equ, Max_deg_phi_direction )
    
    Bp = 20000*10^(-9);     % Equatorial field, T
    Rp = 60280*10^3;        % Planet radius, m
    
    mu_0 = 4*pi*10^(-7);
    Pswnpa = 0.02;
    
    b0 = sqrt(2*mu_0*Pswnpa*10^(-9));   % Field scale
    r0 = (2*Bp*Rp^3/b0)^(1/3);          % Distance scale
    M = (Bp/b0)*(Rp/r0)^3;

    
    
    Nb_points_theta = Max_deg_equ/delta_theta_deg+1;    % Number of points on each direction of the grid, from 0 to Max_deg_equ
    N_theta = Nb_points_theta-2; % Number of points on theta-direction of the grid, from 1 to Max_deg_equ-1

    Nb_points_phi= Max_deg_phi_direction/delta_phi_deg+1;    % Number of points on each direction of the grid, from 0 to Max_deg_equ
    N_phi = Nb_points_phi-2; % Number of points on theta-direction of the grid, from 1 to Max_deg_equ-1

    
    
    Max_deg_theta = Max_deg_equ-delta_theta_deg;
    Max_deg_phi= Max_deg_phi_direction-delta_phi_deg;
    
    delta_theta = delta_theta_deg*(pi/180);
    delta_phi = delta_phi_deg*(pi/180);
    
    theta_vec_1 = repmat(delta_theta_deg:delta_theta_deg:Max_deg_theta,[N_phi 1]);
    theta_vec = theta_vec_1(:)*pi/180;

    phi_vec_1 = repmat(delta_phi_deg:delta_phi_deg:Max_deg_phi,[1 N_theta]).';
    phi_vec = phi_vec_1(:)*pi/180;
    
    
    if k == 1
        
        n_vec = [ 1, -1/r(k)*(r(k)-rss) / (delta_theta), -1/(r(k)*sin(theta_vec(k))) * (r(k+1)-r_eq(2)) / (2*delta_phi)];
        B = M*(1/r(k))^3 * [ -2*sin(theta_vec(k))*sin(phi_vec(k)), cos(theta_vec(k))*sin(phi_vec(k)), cos(phi_vec(k)) ];
        v = [ -cos(theta_vec(k)), sin(theta_vec(k)), 0 ];        
        
    elseif k == N_phi
        
       n_vec = [ 1, -1/r(k)*(r(k)-rss) / (delta_theta), -1/(r(k)*sin(theta_vec(k))) * (r_mer(2)-r(k-1)) / (2*delta_phi)];
        B = M*(1/r(k))^3 * [ -2*sin(theta_vec(k))*sin(phi_vec(k)), cos(theta_vec(k))*sin(phi_vec(k)), cos(phi_vec(k)) ];
        v = [ -cos(theta_vec(k)), sin(theta_vec(k)), 0 ];
        
    elseif k == N_phi*N_theta-N_phi+1

        n_vec = [ 1, -1/r(k)*(r(k) -r(k-N_phi))/ (delta_theta), -1/(r(k)*sin(theta_vec(k))) * (r(k+1)-r_eq(N_theta+1)) / (2*delta_phi)];
        B = M*(1/r(k))^3 * [ -2*sin(theta_vec(k))*sin(phi_vec(k)), cos(theta_vec(k))*sin(phi_vec(k)), cos(phi_vec(k)) ];
        v = [ -cos(theta_vec(k)), sin(theta_vec(k)), 0 ];        
        
    elseif k == N_phi*N_theta
        
         n_vec = [ 1, -1/r(k)*(r(k)-r(k-N_phi))/ (delta_theta), -1/(r(k)*sin(theta_vec(k))) * (r_mer(end-1)-r(k-1)) / (2*delta_phi)];
        B = M*(1/r(k))^3 * [ -2*sin(theta_vec(k))*sin(phi_vec(k)), cos(theta_vec(k))*sin(phi_vec(k)), cos(phi_vec(k)) ];
        v = [ -cos(theta_vec(k)), sin(theta_vec(k)), 0 ];       
        
    elseif k>1 && k<N_phi   %%% left boundary

        n_vec = [ 1, -1/r(k)*(r(k)-rss) / (delta_theta), -1/(r(k)*sin(theta_vec(k))) * (r(k+1)-r(k-1)) / (2*delta_phi)];
        B = M*(1/r(k))^3 * [ -2*sin(theta_vec(k))*sin(phi_vec(k)), cos(theta_vec(k))*sin(phi_vec(k)), cos(phi_vec(k)) ];
        v = [ -cos(theta_vec(k)), sin(theta_vec(k)), 0 ];
        
    elseif k>1 && mod(k,N_phi)==1 && k<N_phi*N_theta-N_phi+1    %%% bottom boundary

        n_vec = [ 1, -1/r(k)*(r(k)-r(k-N_phi)) / (delta_theta), -1/(r(k)*sin(theta_vec(k))) * (r(k+1)-r_eq((k-1)/N_phi+2)) / (2*delta_phi)];
        B = M*(1/r(k))^3 * [ -2*sin(theta_vec(k))*sin(phi_vec(k)), cos(theta_vec(k))*sin(phi_vec(k)), cos(phi_vec(k)) ];
        v = [ -cos(theta_vec(k)), sin(theta_vec(k)), 0 ];        
        
    elseif k>N_phi && mod(k,N_phi)==0 && k<N_phi*N_theta  %%% upper boundary
        
        n_vec = [ 1, -1/r(k)*(r(k)-r(k-N_phi)) / (delta_theta), -1/(r(k)*sin(theta_vec(k))) * (r_mer(k/N_phi+1)-r(k-1)) / (2*delta_phi)];
        B = M*(1/r(k))^3 * [ -2*sin(theta_vec(k))*sin(phi_vec(k)), cos(theta_vec(k))*sin(phi_vec(k)), cos(phi_vec(k)) ];
        v = [ -cos(theta_vec(k)), sin(theta_vec(k)), 0 ];        
        
    elseif k>N_phi*N_theta-N_phi+1 && k<N_phi*N_theta  %%% right boundary

        n_vec = [ 1, -1/r(k)*(r(k)-r(k-N_phi)) / (delta_theta), -1/(r(k)*sin(theta_vec(k))) * (r(k+1)-r(k-1)) / (2*delta_phi)];
        B = M*(1/r(k))^3 * [ -2*sin(theta_vec(k))*sin(phi_vec(k)), cos(theta_vec(k))*sin(phi_vec(k)), cos(phi_vec(k)) ];
        v = [ -cos(theta_vec(k)), sin(theta_vec(k)), 0 ];
        
    else

        n_vec = [ 1, -1/r(k)*(r(k)-r(k-N_phi)) / (delta_theta), -1/(r(k)*sin(theta_vec(k))) * (r(k+1)-r(k-1)) / (2*delta_phi)];
        B = M*(1/r(k))^3 * [ -2*sin(theta_vec(k))*sin(phi_vec(k)), cos(theta_vec(k))*sin(phi_vec(k)), cos(phi_vec(k)) ];
        v = [ -cos(theta_vec(k)), sin(theta_vec(k)), 0 ];

    end
    
    PressureBalanceDiscreteFTCS = norm(cross(n_vec, -B)) + (1/2)*dot(n_vec,v);
%     ScalarProduct = (1/2)*dot(n_vec,v);
end