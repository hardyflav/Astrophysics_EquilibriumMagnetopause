function [F_values_Bottom, Jacobian_Bottom] = F_valuesJacobian_Numerical_Bottom(r, r_eq, r_mer, rLeft, Max_deg_equ, Max_deg_phi_direction, delta_theta_deg, delta_phi_deg, Nb_points_grid, SystemParameters, epsilon)

    beta = SystemParameters.beta;
    F_values_Bottom = ones(Nb_points_grid, 1);
    
    for k = 1:Nb_points_grid

          F_values_Bottom(k) = PB_Cartesian_Bottom( k, r, rLeft, r_eq, r_mer, delta_theta_deg, delta_phi_deg, Max_deg_equ, Max_deg_phi_direction, beta, SystemParameters );
          
    end

    Jacobian_Bottom =  NumericalJacobian_Bottom(Max_deg_equ, Max_deg_phi_direction, delta_theta_deg, delta_phi_deg, r, rLeft, r_eq, r_mer, SystemParameters, epsilon);

end
