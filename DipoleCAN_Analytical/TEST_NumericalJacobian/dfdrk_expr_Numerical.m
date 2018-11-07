function dfdrk_expr_Numerical = dfdrk_expr_Numerical(k, rk, rkp1, rkm1, rkpN, rkmN, rss, r_eq, r_mer, delta_theta_deg, delta_phi_deg, Max_deg_equ, Max_deg_phi_direction, SystemParameters, epsilon)

%     dfdrk_expr_Numerical =  ( PB_Cartesian_Numerical( k, rk+epsilon, rkp1, rkm1, rkpN, rkmN, rss, r_eq, r_mer, delta_theta_deg, delta_phi_deg, Max_deg_equ, Max_deg_phi_direction, SystemParameters ) - PB_Cartesian_Numerical( k, rk, rkp1, rkm1, rkpN, rkmN, rss, r_eq, r_mer, delta_theta_deg, delta_phi_deg, Max_deg_equ, Max_deg_phi_direction, SystemParameters ) ) / epsilon;
    dfdrk_expr_Numerical =  ( PB_Cartesian_Numerical( k, rk+epsilon/2, rkp1, rkm1, rkpN, rkmN, rss, r_eq, r_mer, delta_theta_deg, delta_phi_deg, Max_deg_equ, Max_deg_phi_direction, SystemParameters ) - PB_Cartesian_Numerical( k, rk-epsilon/2, rkp1, rkm1, rkpN, rkmN, rss, r_eq, r_mer, delta_theta_deg, delta_phi_deg, Max_deg_equ, Max_deg_phi_direction, SystemParameters ) ) / (epsilon);

end