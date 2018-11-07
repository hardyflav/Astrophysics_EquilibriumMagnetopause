function FValuesFringingField = FValuesFringingField(r, r_eq, r_right_outer_boundary, r_mer, rss, Max_deg_equ, Max_deg_phi_direction, delta_theta_deg, delta_phi_deg, Nb_points_grid, BfPrevious, BIMF)
    
    FValuesFringingField = ones(Nb_points_grid, 1);
    
    for k = 1:Nb_points_grid
        FValuesFringingField(k) = PressureBalanceFringingField( k, r, rss, r_eq, r_right_outer_boundary, r_mer, delta_theta_deg, delta_phi_deg, Max_deg_equ, Max_deg_phi_direction, BfPrevious, BIMF );
    end

end
