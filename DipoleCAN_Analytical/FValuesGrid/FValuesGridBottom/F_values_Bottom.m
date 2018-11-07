function [F_values] = F_values_Bottom(r, r_eq, r_mer, rLeft, ThetaLeftBoundaryDeg, Max_deg_equ, Max_deg_phi_direction, delta_theta_deg, delta_phi_deg, Nb_points_grid, SystemParameters)
    
    %N = Max_deg_theta/delta_theta_deg+1;
    %Nb_points_grid = (N-(1/delta_theta_deg+1)+1)^2;
    
%     
%     [~, ~,...
%     N_theta, N_phi,...
%     ~, ~] = GridDetailsBottom(ThetaLeftBoundaryDeg, Max_deg_equ, Max_deg_phi_direction, delta_theta_deg, delta_phi_deg);

%     rGrid = reshape(r, [N_phi, N_theta]);
%     r_eq_temp = 2*rGrid(2,:)-rGrid(1,:);
%     r_mer_temp = 2*rGrid(end,:)-rGrid(end-1,:);
%     
%     r_eq = horzcat(0, r_eq_temp, 0).';
%     r_mer = horzcat(0, r_mer_temp, 0).';
%     
    
    
    F_values = ones(Nb_points_grid, 1);
    
%     ScalarProducts = ones(Nb_points_grid, 1);
    
    for k = 1:Nb_points_grid
        F_values(k) = PB_Cartesian_Bottom( k, r, rLeft, r_eq, r_mer, delta_theta_deg, delta_phi_deg, ThetaLeftBoundaryDeg, Max_deg_equ, Max_deg_phi_direction, SystemParameters);
    end

end
