function Jacobian = Jacobian(Max_deg_equ, Max_deg_phi_direction, delta_theta_deg, delta_phi_deg, r, rss, r_eq, r_mer)

% r_guess is such that theta,phi = delta:delta:Max-delta, phi
% ie it should contain (N_theta-2)*(N_phi-2) components

    [~, ~,...
    N_theta, N_phi,...
    theta_vec, phi_vec] = GridDetails(Max_deg_equ, Max_deg_phi_direction, delta_theta_deg, delta_phi_deg);

    delta_theta = delta_theta_deg*(pi/180);
    delta_phi = delta_phi_deg*(pi/180);

I=zeros((5*N_phi*N_theta-2*N_phi-2), 1);    % Vector of row indices in the Jacobian matrix J
J=zeros((5*N_phi*N_theta-2*N_phi-2), 1);      % Vector of column indices in the Jacobian matrix J
V=zeros((5*N_phi*N_theta-2*N_phi-2), 1);     % Vector of values in the Jacobian matrix J


iter=1;
for i=1:N_phi*N_theta
    
    for j=1:N_phi*N_theta
        
        if i==j            
           I(iter)=i;
           J(iter)=j;
           V(iter)=dfkdrk(i,Max_deg_equ, Max_deg_phi_direction, theta_vec, phi_vec, delta_theta, delta_phi, r, rss, r_eq, r_mer);
           iter=iter+1;
        
        
        elseif j==i+1       % sup-diagonal
            I(iter)=i;
            J(iter)=j;            
            V(iter)=dfkdrkp1(i,Max_deg_equ, Max_deg_phi_direction, theta_vec, phi_vec, delta_theta, delta_phi, r, rss, r_eq);
            iter=iter+1;
            
        elseif j==i-1  % sub-diagonal
            I(iter)=i;
            J(iter)=j;            
            V(iter)= dfkdrkm1(i,Max_deg_equ, Max_deg_phi_direction, theta_vec, phi_vec, delta_theta, delta_phi, r, rss, r_mer);
            iter=iter+1;            
        
        elseif j==i-N_phi  % N-sub-diagonal
            I(iter)=i;
            J(iter)=j;            
            V(iter)= dfkdrkmN(i,Max_deg_equ, Max_deg_phi_direction, theta_vec, phi_vec, delta_theta, delta_phi, r, r_eq, r_mer);
            iter=iter+1;     

        elseif j==i+N_phi      % N-sup-diagonal
           I(iter)=i;
           J(iter)=j;            
           V(iter)=dfkdrkpN(i,Max_deg_equ, Max_deg_phi_direction, theta_vec, phi_vec, delta_theta, delta_phi, r, rss, r_eq, r_mer);
           iter=iter+1;
        
        end
        
    end    
end

Jacobian = sparse(I,J,V);

end