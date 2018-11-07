function [F_vec, J] = Values_Jacobian( r, r_eq, rss, r_mer, Max_deg_equ, Max_deg_phi, delta_theta_deg, delta_phi_deg, Nb_points_grid)

    F_vec = F_values(r, r_eq, r_mer, rss, Max_deg_equ, Max_deg_phi, delta_theta_deg, delta_phi_deg, Nb_points_grid);
    J = Jacobian(Max_deg_equ, Max_deg_phi, delta_theta_deg, delta_phi_deg, r, rss, r_eq, r_mer);

end