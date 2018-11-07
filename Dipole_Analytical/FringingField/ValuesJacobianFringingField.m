function [F_vec, J] = ValuesJacobianFringingField( r, r_eq, rss, r_right_outer_boundary, r_mer, Max_deg_equ, Max_deg_phi, delta_theta_deg, delta_phi_deg, Nb_points_grid, BfPrevious, BIMF)

    F_vec = FValuesFringingField(r, r_eq, r_right_outer_boundary, r_mer, rss, Max_deg_equ, Max_deg_phi, delta_theta_deg, delta_phi_deg, Nb_points_grid, BfPrevious, BIMF);
    J = Jacobian(Max_deg_equ, Max_deg_phi, delta_theta_deg, delta_phi_deg, r, rss, r_eq, r_right_outer_boundary, r_mer);

end