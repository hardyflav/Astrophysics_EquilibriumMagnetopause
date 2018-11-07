function [F_vec, J] = Values_Jacobian_Bottom( r, r_eq, rLeft, r_right_outer_boundary, r_mer, ThetaLeftBoundaryDeg, Max_deg_equ, Max_deg_phi, delta_theta_deg, delta_phi_deg, Nb_points_grid)

    F_vec = F_values_Bottom(r, r_eq, r_right_outer_boundary, r_mer, rLeft, ThetaLeftBoundaryDeg, Max_deg_equ, Max_deg_phi, delta_theta_deg, delta_phi_deg, Nb_points_grid);
    J = Jacobian_Bottom(ThetaLeftBoundaryDeg, Max_deg_equ, Max_deg_phi, delta_theta_deg, delta_phi_deg, r, rLeft, r_eq, r_right_outer_boundary, r_mer);

end