function JacobianPattern = JacobianPattern(Max_deg_equ, Max_deg_phi_direction, delta_theta_deg, delta_phi_deg)

% r_guess is such that theta,phi = delta:delta:Max-delta, phi
% ie it should contain (N_theta-2)*(N_phi-2) components


    Nb_points_theta = Max_deg_equ/delta_theta_deg+1;    % Number of points on each direction of the grid, from 0 to Max_deg_equ
    N_theta = Nb_points_theta-2; % Number of points on theta-direction of the grid, from 1 to Max_deg_equ-1

    Nb_points_phi= Max_deg_phi_direction/delta_phi_deg+1;    % Number of points on each direction of the grid, from 0 to Max_deg_equ
    N_phi = Nb_points_phi-2; % Number of points on theta-direction of the grid, from 1 to Max_deg_equ-1

    
    
I=zeros(N_phi*N_theta, 1);      % Vector of row indices in the Jacobian matrix J
J=zeros(N_phi*N_theta, 1);      % Vector of column indices in the Jacobian matrix J
V=zeros(N_phi*N_theta, 1);     % Vector of values in the Jacobian matrix J


iter=1;
for i=1:N_phi*N_theta
    
    for j=1:N_phi*N_theta
        
        if i==j            
           I(iter)=i;
           J(iter)=j;
           V(iter)=1;
           iter=iter+1;
        
        
        elseif j==i+1       % sup-diagonal
            I(iter)=i;
            J(iter)=j;            
            V(iter)=1;
            iter=iter+1;
            
        elseif j==i-1  % sub-diagonal
            I(iter)=i;
            J(iter)=j;            
            V(iter)= 1;
            iter=iter+1;            
        
        elseif j==i-N_phi  % N-sub-diagonal
            I(iter)=i;
            J(iter)=j;            
            V(iter)= 1;
            iter=iter+1;   
     

        elseif j==i+N_phi      % N-sup-diagonal
           I(iter)=i;
           J(iter)=j;            
           V(iter)=1;
           iter=iter+1;
           
           
        %else
            %I(iter)=i;
            %J(iter)=j;
            %iter=iter+1;
        
        end
        
    end    
end

JacobianPattern = sparse(I,J,V);

end