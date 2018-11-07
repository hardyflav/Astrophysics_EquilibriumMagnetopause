
function dfkdrkpNBottom = dfkdrkpN_Bottom(k, ThetaLeftBoundaryDeg, Max_deg_equ, Max_deg_phi_direction, theta_vec, phi_vec, delta_theta, delta_phi, r, rLeft, r_eq, r_mer)

  [~, ~,...
    N_theta, N_phi,...
    ~, ~] = GridDetailsBottom(ThetaLeftBoundaryDeg, Max_deg_equ, Max_deg_phi_direction, delta_theta*180/pi, delta_phi*180/pi);

    rLeftBoundary = rLeft(2:end-1);


    if k==1
        
        dfkdrkpNBottom = dfdrkpN_expr(k,theta_vec, phi_vec, delta_theta, delta_phi, r(k), r(k+1), r_eq( (k-1)/N_phi+2 ), r(k+N_phi), rLeftBoundary(1));
    
    elseif k==N_phi*N_theta
        
        dfkdrkpNBottom = 0;
        
    elseif k==N_phi
        
        dfkdrkpNBottom = dfdrkpN_expr(k,theta_vec, phi_vec, delta_theta, delta_phi, r(k), r_mer(2), r(k-1), r(k+N_phi), rLeftBoundary(end));
                 
    elseif k==N_phi*N_theta-N_phi+1
        
        dfkdrkpNBottom = 0;
        
    elseif k>N_phi*N_theta-N_phi+1 && k<N_phi*N_theta  %%% right boundary
        
        dfkdrkpNBottom = 0;
        
    elseif k>1 && mod(k,N_phi)==1 && k<N_phi*N_theta-N_phi+1    %%% bottom boundary
        
        dfkdrkpNBottom = dfdrkpN_expr(k,theta_vec, phi_vec, delta_theta, delta_phi, r(k), r(k+1), r_eq( (k-1)/N_phi+2 ), r(k+N_phi), r(k-N_phi));
         
    elseif k>N_phi && mod(k,N_phi)==0 && k<N_phi*N_theta  %%% upper boundary
        
            dfkdrkpNBottom = dfdrkpN_expr(k,theta_vec, phi_vec, delta_theta, delta_phi, r(k), r_mer(k/N_phi+1), r(k-1), r(k+N_phi), r(k-N_phi));
        
    elseif k>1 && k<N_phi   %%% left boundary
        
        dfkdrkpNBottom = dfdrkpN_expr(k,theta_vec, phi_vec, delta_theta, delta_phi, r(k), r(k+1), r(k-1), r(k+N_phi), rLeftBoundary(k));
                          
    else
        
        dfkdrkpNBottom = dfdrkpN_expr(k,theta_vec, phi_vec, delta_theta, delta_phi, r(k), r(k+1), r(k-1), r(k+N_phi), r(k-N_phi));
        
    end
end