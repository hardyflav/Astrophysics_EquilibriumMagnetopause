
function dfkdrkm1Bottom = dfkdrkm1_Bottom(k, ThetaLeftBoundaryDeg, Max_deg_equ, Max_deg_phi_direction, theta_vec, phi_vec, delta_theta, delta_phi, r, rLeft, r_eq, r_mer)

  
  [~, ~,...
    N_theta, N_phi,...
    ~, ~] = GridDetailsBottom(ThetaLeftBoundaryDeg, Max_deg_equ, Max_deg_phi_direction, delta_theta*180/pi, delta_phi*180/pi);

    rLeftBoundary = rLeft(2:end-1);


    if k==1
        
        dfkdrkm1Bottom = 0;
        
    elseif k==N_phi*N_theta
        
        dfkdrkm1Bottom = dfkdrkm1Right_expr(k ,theta_vec, phi_vec, delta_theta, delta_phi, r(k), r_mer(N_theta+1), r(k-1), r(k-N_phi));
        
    elseif k==N_phi
        
        dfkdrkm1Bottom = dfdrkm1_expr(k, theta_vec, phi_vec, delta_theta, delta_phi, r(k), r_mer(2), r(k-1), r(k+N_phi), rLeftBoundary(end));
      
    elseif k==N_phi*N_theta-N_phi+1
        
        dfkdrkm1Bottom = 0;
        
    elseif k>N_phi*N_theta-N_phi+1 && k<N_phi*N_theta  %%% right boundary
              
        dfkdrkm1Bottom = dfkdrkm1Right_expr(k, theta_vec, phi_vec, delta_theta, delta_phi, r(k), r(k+1), r(k-1), r(k-N_phi));
                
    elseif k>1 && mod(k,N_phi)==1 && k<N_phi*N_theta-N_phi+1    %%% bottom boundary

        dfkdrkm1Bottom = 0;
        
    elseif k>N_phi && mod(k,N_phi)==0 && k<N_phi*N_theta  %%% upper boundary

        dfkdrkm1Bottom = dfdrkm1_expr(k, theta_vec, phi_vec, delta_theta, delta_phi, r(k), r_mer(k/N_phi+1), r(k-1), r(k+N_phi), r(k-N_phi));
            
    elseif k>1 && k<N_phi   %%% left boundary     
        
        dfkdrkm1Bottom = dfdrkm1_expr(k, theta_vec, phi_vec, delta_theta, delta_phi, r(k), r(k+1), r(k-1), r(k+N_phi), rLeftBoundary(k));        
         
    else
        
        dfkdrkm1Bottom = dfdrkm1_expr(k, theta_vec, phi_vec, delta_theta, delta_phi, r(k), r(k+1), r(k-1), r(k+N_phi), r(k-N_phi));
        
    end
    
 end