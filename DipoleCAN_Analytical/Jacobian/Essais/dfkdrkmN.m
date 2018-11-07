
function dfkdrkmN = dfkdrkmN(k,Max_deg_equ, Max_deg_phi_direction, theta_vec, phi_vec, delta_theta, delta_phi, r, r_eq, r_mer)
    
    [~, ~,...
    N_theta, N_phi,...
    ~, ~] = GridDetails(Max_deg_equ, Max_deg_phi_direction, delta_theta*180/pi, delta_phi*180/pi);


    if k==1
        
        dfkdrkmN = 0;    
    
    elseif k==N_phi*N_theta
        
        dfkdrkmN = dfkdrkmNRight_expr(k,theta_vec, phi_vec, delta_theta, delta_phi, r(k), r_mer(N_theta+1), r(k-1), r(k-N_phi));
        
    elseif k==N_phi
        
        dfkdrkmN = 0;
         
    elseif k==N_phi*N_theta-N_phi+1
        
        dfkdrkmN = dfkdrkmNRight_expr(k,theta_vec, phi_vec, delta_theta, delta_phi, r(k), r(k+1), r_eq( N_theta+1 ), r(k-N_phi));
         
    elseif k>N_phi*N_theta-N_phi+1 && k<N_phi*N_theta  %%% right boundary
        
        dfkdrkmN = dfkdrkmNRight_expr(k,theta_vec, phi_vec, delta_theta, delta_phi, r(k), r(k+1), r(k-1), r(k-N_phi));
        
    elseif k>1 && mod(k,N_phi)==1 && k<N_phi*N_theta-N_phi+1    %%% bottom boundary
        
        dfkdrkmN = dfdrkmN_expr(k,theta_vec, phi_vec, delta_theta, delta_phi, r(k), r(k+1), r_eq( (k-1)/N_phi+2 ), r(k+N_phi), r(k-N_phi));
         
    elseif k>N_phi && mod(k,N_phi)==0 && k<N_phi*N_theta  %%% upper boundary
        
       dfkdrkmN = dfdrkmN_expr(k,theta_vec, phi_vec, delta_theta, delta_phi, r(k), r_mer(k/N_phi+1), r(k-1), r(k+N_phi), r(k-N_phi));
 
    elseif k>1 && k<N_phi   %%% left boundary
        
       dfkdrkmN = 0;          
        
    else
        
        dfkdrkmN = dfdrkmN_expr(k, theta_vec, phi_vec, delta_theta, delta_phi, r(k), r(k+1), r(k-1), r(k+N_phi), r(k-N_phi));
       
    end
end