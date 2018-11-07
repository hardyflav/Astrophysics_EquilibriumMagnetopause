function [F_values] = F_values(r, r_eq, r_mer, rss, Max_deg_equ, Max_deg_phi_direction, delta_theta_deg, delta_phi_deg, Nb_points_grid, beta, SystemParameters)

    F_values = ones(Nb_points_grid, 1);
    
    for k = 1:Nb_points_grid

%         F_values(k) = pressure_balance_v3( k, r, rss, r_eq, r_mer, delta_theta_deg, delta_phi_deg, Max_deg_equ, Max_deg_phi_direction, beta, CAN_DiskParameters );
          F_values(k) = PB_Cartesian( k, r, rss, r_eq, r_mer, delta_theta_deg, delta_phi_deg, Max_deg_equ, Max_deg_phi_direction, beta, SystemParameters );
          
    end

end
