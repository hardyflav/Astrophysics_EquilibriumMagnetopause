
function dfkdrkmN_Numerical_Bottom = dfkdrkmN_Numerical_Bottom(k,Max_deg_equ, Max_deg_phi_direction, delta_theta, delta_phi, r, rLeft, r_eq, r_mer, ThetaLeftBoundaryDeg, SystemParameters, epsilon)

%     Nb_points_theta = size(r_eq, 1);    % Number of points on each direction of the grid, from ThetaLeftBoundaryDeg to Max_deg_equ
%     N_theta = Nb_points_theta-2; % Number of points on theta-direction of the grid, from 1 to Max_deg_equ-1
% 
%     Nb_points_phi= size(rLeft, 1);    % Number of points on each direction of the grid, from 0 to Max_deg_equ
%     N_phi = Nb_points_phi-2; % Number of points on theta-direction of the grid, from 1 to Max_deg_equ-1
    
    delta_theta_deg = delta_theta*180/pi;
    delta_phi_deg = delta_phi*180/pi;

    [~, ~,...
    N_theta, N_phi,...
    ~, ~] = GridDetailsBottom(ThetaLeftBoundaryDeg, Max_deg_equ, Max_deg_phi_direction, delta_theta_deg, delta_phi_deg);
    
    rLeftBoundary = rLeft(2:end-1);
    
    if k==1

        dfkdrkmN_Numerical_Bottom = dfdrkmN_expr_Numerical_Bottom(k,  r(k), r(k+1), r_eq(2), r(k+N_phi), rLeftBoundary(k), rLeft, r_eq, r_mer, delta_theta_deg, delta_phi_deg, Max_deg_equ, Max_deg_phi_direction, ThetaLeftBoundaryDeg, SystemParameters, epsilon);
         
    elseif k==N_phi*N_theta
    
        dfkdrkmN_Numerical_Bottom = dfdrkmN_expr_Numerical_Bottom(k,  r(k), r_mer(N_theta+1), r(k-1), 0, r(k-N_phi), rLeft, r_eq, r_mer, delta_theta_deg, delta_phi_deg, Max_deg_equ, Max_deg_phi_direction, ThetaLeftBoundaryDeg, SystemParameters, epsilon);
        
    elseif k==N_phi
        
        dfkdrkmN_Numerical_Bottom = dfdrkmN_expr_Numerical_Bottom(k,  r(k), r_mer(2), r(k-1), r(k+N_phi), rLeftBoundary(k), rLeft, r_eq, r_mer, delta_theta_deg, delta_phi_deg, Max_deg_equ, Max_deg_phi_direction, ThetaLeftBoundaryDeg, SystemParameters, epsilon);
           
    elseif k==N_phi*N_theta-N_phi+1        

        dfkdrkmN_Numerical_Bottom = dfdrkmN_expr_Numerical_Bottom(k,  r(k), r(k+1), r_eq(N_theta+1), 0, r(k-N_phi), rLeft, r_eq, r_mer, delta_theta_deg, delta_phi_deg, Max_deg_equ, Max_deg_phi_direction, ThetaLeftBoundaryDeg, SystemParameters, epsilon);
        
    elseif k>N_phi*N_theta-N_phi+1 && k<N_phi*N_theta  %%% right boundary
        
        dfkdrkmN_Numerical_Bottom = dfdrkmN_expr_Numerical_Bottom(k, r(k), r(k+1), r(k-1), 0, r(k-N_phi), rLeft, r_eq, r_mer, delta_theta_deg, delta_phi_deg, Max_deg_equ, Max_deg_phi_direction, ThetaLeftBoundaryDeg, SystemParameters, epsilon);
        
    elseif k>1 && mod(k,N_phi)==1 && k<N_phi*N_theta-N_phi+1    %%% bottom boundary

        dfkdrkmN_Numerical_Bottom = dfdrkmN_expr_Numerical_Bottom(k, r(k), r(k+1), r_eq( (k-1)/N_phi+2 ), r(k+N_phi), r(k-N_phi), rLeft, r_eq, r_mer, delta_theta_deg, delta_phi_deg, Max_deg_equ, Max_deg_phi_direction, ThetaLeftBoundaryDeg, SystemParameters, epsilon);
        
    elseif k>N_phi && mod(k,N_phi)==0 && k<N_phi*N_theta  %%% upper boundary
 
        dfkdrkmN_Numerical_Bottom = dfdrkmN_expr_Numerical_Bottom(k, r(k), r_mer(k/N_phi+1), r(k-1), r(k+N_phi), r(k-N_phi), rLeft, r_eq, r_mer, delta_theta_deg, delta_phi_deg, Max_deg_equ, Max_deg_phi_direction, ThetaLeftBoundaryDeg, SystemParameters, epsilon);
            
    elseif k>1 && k<N_phi   %%% left boundary
        
        dfkdrkmN_Numerical_Bottom = dfdrkmN_expr_Numerical_Bottom(k, r(k), r(k+1), r(k-1), r(k+N_phi), rLeftBoundary(k), rLeft, r_eq, r_mer, delta_theta_deg, delta_phi_deg, Max_deg_equ, Max_deg_phi_direction, ThetaLeftBoundaryDeg, SystemParameters, epsilon);
     
    else
        
        dfkdrkmN_Numerical_Bottom = dfdrkmN_expr_Numerical_Bottom(k, r(k), r(k+1), r(k-1), r(k+N_phi), r(k-N_phi), rLeft, r_eq, r_mer, delta_theta_deg, delta_phi_deg, Max_deg_equ, Max_deg_phi_direction, ThetaLeftBoundaryDeg, SystemParameters, epsilon);
         
    end
end

