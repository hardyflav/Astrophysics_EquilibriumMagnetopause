function [pressure_balance_v3] = pressure_balance_v3( k, r, rss, r_eq, r_mer, delta_theta_deg, delta_phi_deg, Max_deg_equ, Max_deg_phi_direction )
    
    Bp = 20000*10^(-9);     % Equatorial field, T
    Rp = 60280*10^3;        % Planet radius, m
    
    mu_0 = 4*pi*10^(-7);
    Pswnpa = 0.02;
    
    b0 = sqrt(2*mu_0*Pswnpa*10^(-9));   % Field scale
    r0 = (2*Bp*Rp^3/b0)^(1/3);          % Distance scale
    M = (Bp/b0)*(Rp/r0)^3;

    [~, ~,...
    N_theta, N_phi,...
    theta_vec, phi_vec] = GridDetails(Max_deg_equ, Max_deg_phi_direction, delta_theta_deg, delta_phi_deg);

    delta_theta = delta_theta_deg*(pi/180);
    delta_phi = delta_phi_deg*(pi/180);
    
    if k == 1
        
        n_vec = [ 1, -1/r(k)*(r(k+N_phi)-rss) / (2*delta_theta), -1/(r(k)*sin(theta_vec(k))) * (r(k+1)-r_eq(2)) / (2*delta_phi)];
        B = M*(1/r(k))^3 * [ -2*sin(theta_vec(k))*sin(phi_vec(k)), cos(theta_vec(k))*sin(phi_vec(k)), cos(phi_vec(k)) ];
        v = [ -cos(theta_vec(k)), sin(theta_vec(k)), 0 ];        
        
    elseif k == N_phi
        
        n_vec = [ 1, -1/r(k)*(r(k+N_phi)-rss) / (2*delta_theta), -1/(r(k)*sin(theta_vec(k))) * (r_mer(2)-r(k-1)) / (2*delta_phi)];
        B = M*(1/r(k))^3 * [ -2*sin(theta_vec(k))*sin(phi_vec(k)), cos(theta_vec(k))*sin(phi_vec(k)), cos(phi_vec(k)) ];
        v = [ -cos(theta_vec(k)), sin(theta_vec(k)), 0 ];
        
    elseif k == N_phi*N_theta-N_phi+1

        n_vec = [ 1, -1/r(k)*(2*r(k)-r(k-N_phi)-r(k-N_phi))/ (2*delta_theta), -1/(r(k)*sin(theta_vec(k))) * (r(k+1)-r_eq(N_theta+1)) / (2*delta_phi)];
        B = M*(1/r(k))^3 * [ -2*sin(theta_vec(k))*sin(phi_vec(k)), cos(theta_vec(k))*sin(phi_vec(k)), cos(phi_vec(k)) ];
        v = [ -cos(theta_vec(k)), sin(theta_vec(k)), 0 ];        
        
    elseif k == N_phi*N_theta
        
        n_vec = [ 1, -1/r(k)*(2*r(k)-r(k-N_phi)-r(k-N_phi))/ (2*delta_theta), -1/(r(k)*sin(theta_vec(k))) * (r_mer(N_theta+1)-r(k-1)) / (2*delta_phi)];
        B = M*(1/r(k))^3 * [ -2*sin(theta_vec(k))*sin(phi_vec(k)), cos(theta_vec(k))*sin(phi_vec(k)), cos(phi_vec(k)) ];
        v = [ -cos(theta_vec(k)), sin(theta_vec(k)), 0 ];       
        
    elseif k>1 && k<N_phi   %%% left boundary

        n_vec = [ 1, -1/r(k)*(r(k+N_phi)-rss) / (2*delta_theta), -1/(r(k)*sin(theta_vec(k))) * (r(k+1)-r(k-1)) / (2*delta_phi)];
        B = M*(1/r(k))^3 * [ -2*sin(theta_vec(k))*sin(phi_vec(k)), cos(theta_vec(k))*sin(phi_vec(k)), cos(phi_vec(k)) ];
        v = [ -cos(theta_vec(k)), sin(theta_vec(k)), 0 ];
        
    elseif k>1 && mod(k,N_phi)==1 && k<N_phi*N_theta-N_phi+1    %%% bottom boundary

        n_vec = [ 1, -1/r(k)*(r(k+N_phi)-r(k-N_phi)) / (2*delta_theta), -1/(r(k)*sin(theta_vec(k))) * (r(k+1)-r_eq((k-1)/N_phi+2)) / (2*delta_phi)];
        B = M*(1/r(k))^3 * [ -2*sin(theta_vec(k))*sin(phi_vec(k)), cos(theta_vec(k))*sin(phi_vec(k)), cos(phi_vec(k)) ];
        v = [ -cos(theta_vec(k)), sin(theta_vec(k)), 0 ];        
        
    elseif k>N_phi && mod(k,N_phi)==0 && k<N_phi*N_theta  %%% upper boundary
        
        n_vec = [ 1, -1/r(k)*(r(k+N_phi)-r(k-N_phi)) / (2*delta_theta), -1/(r(k)*sin(theta_vec(k))) * (r_mer(k/N_phi+1)-r(k-1)) / (2*delta_phi)];
        B = M*(1/r(k))^3 * [ -2*sin(theta_vec(k))*sin(phi_vec(k)), cos(theta_vec(k))*sin(phi_vec(k)), cos(phi_vec(k)) ];
        v = [ -cos(theta_vec(k)), sin(theta_vec(k)), 0 ];        
        
    elseif k>N_phi*N_theta-N_phi+1 && k<N_phi*N_theta  %%% right boundary

        n_vec = [ 1, -1/r(k)*(2*r(k)-r(k-N_phi)-r(k-N_phi)) / (2*delta_theta), -1/(r(k)*sin(theta_vec(k))) * (r(k+1)-r(k-1)) / (2*delta_phi)];
        B = M*(1/r(k))^3 * [ -2*sin(theta_vec(k))*sin(phi_vec(k)), cos(theta_vec(k))*sin(phi_vec(k)), cos(phi_vec(k)) ];
        v = [ -cos(theta_vec(k)), sin(theta_vec(k)), 0 ];
        
    else

        n_vec = [ 1, -1/r(k)*(r(k+N_phi)-r(k-N_phi)) / (2*delta_theta), -1/(r(k)*sin(theta_vec(k))) * (r(k+1)-r(k-1)) / (2*delta_phi)];
        B = M*(1/r(k))^3 * [ -2*sin(theta_vec(k))*sin(phi_vec(k)), cos(theta_vec(k))*sin(phi_vec(k)), cos(phi_vec(k)) ];
        v = [ -cos(theta_vec(k)), sin(theta_vec(k)), 0 ];

    end
    
    pressure_balance_v3 = norm(cross(n_vec, -B)) + (1/2)*dot(n_vec,v);
    
end
