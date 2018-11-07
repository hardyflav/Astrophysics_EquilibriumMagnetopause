
function dfkdrkp1 = dfkdrkp1(k,Max_deg_equ, Max_deg_phi_direction, theta_vec, phi_vec, delta_theta, delta_phi, r, rss, r_eq)
  
     [~, ~,...
    N_theta, N_phi,...
    ~, ~] = GridDetails(Max_deg_equ, Max_deg_phi_direction, delta_theta*180/pi, delta_phi*180/pi);

    if k==1
        
        dfkdrkp1 = dfdrkp1_expr(k,theta_vec, phi_vec, delta_theta, delta_phi, r(k), r(k+1), r_eq( (k-1)/N_phi+2 ), r(k+N_phi), rss);
        
    elseif k==N_phi*N_theta
        
        dfkdrkp1 = 0;
        
    elseif k==N_phi
        
        dfkdrkp1 = 0;
%          
    elseif k==N_phi*N_theta-N_phi+1
        
        dfkdrkp1 = dfkdrkp1Right_expr(k,theta_vec, phi_vec, delta_theta, delta_phi, r(k), r(k+1), r_eq( N_theta+1 ),  r(k-N_phi));
        
    elseif k>N_phi*N_theta-N_phi+1 && k<N_phi*N_theta  %%% right boundary
        
        dfkdrkp1 = dfkdrkp1Right_expr(k,theta_vec, phi_vec, delta_theta, delta_phi, r(k), r(k+1), r(k-1), r(k-N_phi));
        
    elseif k>1 && mod(k,N_phi)==1 && k<N_phi*N_theta-N_phi+1    %%% bottom boundary

        dfkdrkp1 = dfdrkp1_expr(k,theta_vec, phi_vec, delta_theta, delta_phi, r(k), r(k+1), r_eq( (k-1)/N_phi+2 ), r(k+N_phi), r(k-N_phi));
        
    elseif k>N_phi && mod(k,N_phi)==0 && k<N_phi*N_theta  %%% upper boundary
 
        dfkdrkp1 = 0;
        
    elseif k>1 && k<N_phi   %%% left boundary
        
        dfkdrkp1 = dfdrkp1_expr(k,theta_vec, phi_vec, delta_theta, delta_phi, r(k), r(k+1), r(k-1), r(k+N_phi), rss);
                 
    else
        
        dfkdrkp1 = dfdrkp1_expr(k,theta_vec, phi_vec, delta_theta, delta_phi, r(k), r(k+1), r(k-1), r(k+N_phi), r(k-N_phi));
       
    end
    
 end