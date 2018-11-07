function dfdrkp1_expr_Numerical = dfdrkp1_expr_Numerical(k, rk, rkp1, rkm1, rkpN, rkmN, rss, r_eq, r_mer, delta_theta_deg, delta_phi_deg, Max_deg_equ, Max_deg_phi_direction, SystemParameters, epsilon)

%     dfdrkp1_expr_Numerical =  ( PB_Cartesian_Numerical( k, rk, rkp1+epsilon, rkm1, rkpN, rkmN, rss, r_eq, r_mer, delta_theta_deg, delta_phi_deg, Max_deg_equ, Max_deg_phi_direction, SystemParameters ) - PB_Cartesian_Numerical( k, rk, rkp1, rkm1, rkpN, rkmN, rss, r_eq, r_mer, delta_theta_deg, delta_phi_deg, Max_deg_equ, Max_deg_phi_direction, SystemParameters ) ) / epsilon;
    dfdrkp1_expr_Numerical =  ( PB_Cartesian_Numerical( k, rk, rkp1+epsilon/2, rkm1, rkpN, rkmN, rss, r_eq, r_mer, delta_theta_deg, delta_phi_deg, Max_deg_equ, Max_deg_phi_direction, SystemParameters ) - PB_Cartesian_Numerical( k, rk, rkp1-epsilon/2, rkm1, rkpN, rkmN, rss, r_eq, r_mer, delta_theta_deg, delta_phi_deg, Max_deg_equ, Max_deg_phi_direction, SystemParameters ) ) / epsilon;

end