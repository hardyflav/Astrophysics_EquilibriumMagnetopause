function F_valuesRelative = F_valuesRelative(r, r_eq, r_mer, rss, Max_deg_equ, Max_deg_phi_direction, delta_theta_deg, delta_phi_deg, Nb_points_grid)
    
    F_valuesRelative = ones(Nb_points_grid, 1);
    
    for k = 1:Nb_points_grid

        F_valuesRelative(k) = pressure_balance_v3Relative( k, r, rss, r_eq, r_mer, delta_theta_deg, delta_phi_deg, Max_deg_equ, Max_deg_phi_direction );

    end

end
