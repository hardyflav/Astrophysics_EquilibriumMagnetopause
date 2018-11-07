function [F_values] = F_values_Bottom(r, r_eq, r_right_outer_boundary, r_mer, rLeft, ThetaLeftBoundaryDeg, Max_deg_equ, Max_deg_phi_direction, delta_theta_deg, delta_phi_deg, Nb_points_grid)
    
    %N = Max_deg_theta/delta_theta_deg+1;
    %Nb_points_grid = (N-(1/delta_theta_deg+1)+1)^2;

    F_values = ones(Nb_points_grid, 1);
%     ScalarProducts = ones(Nb_points_grid, 1);
    
    for k = 1:Nb_points_grid
        F_values(k) = pressure_balance_v3_Bottom( k, r, rLeft, r_eq, r_right_outer_boundary, r_mer, delta_theta_deg, delta_phi_deg, ThetaLeftBoundaryDeg, Max_deg_equ, Max_deg_phi_direction );
%         [ak] = PressureBalanceDiscreteFTCS( k, r, rss, r_eq, r_mer, delta_theta_deg, delta_phi_deg, Max_deg_equ, Max_deg_phi_direction );
%         ScalarProducts(k) = bk;
    end

end
