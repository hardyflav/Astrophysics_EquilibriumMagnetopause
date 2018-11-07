function Jacobian_Numerical_Bottom = Jacobian_Numerical_Bottom(ThetaLeftBoundaryDeg, Max_deg_equ, Max_deg_phi_direction, delta_theta_deg, delta_phi_deg, r, rLeft, r_eq, r_mer, SystemParameters, epsilon)

% r_guess is such that theta,phi = delta:delta:Max-delta, phi
% ie it should contain (N_theta-2)*(N_phi-2) components

    [~, ~,...
    N_theta, N_phi,...
    ~, ~] = GridDetailsBottom(ThetaLeftBoundaryDeg, Max_deg_equ, Max_deg_phi_direction, delta_theta_deg, delta_phi_deg);

    delta_theta = delta_theta_deg*(pi/180);
    delta_phi = delta_phi_deg*(pi/180);

I=zeros((5*N_phi*N_theta-2*N_phi-2), 1);    % Vector of row indices in the Jacobian matrix J
J=zeros((5*N_phi*N_theta-2*N_phi-2), 1);      % Vector of column indices in the Jacobian matrix J
V=zeros((5*N_phi*N_theta-2*N_phi-2), 1);     % Vector of values in the Jacobian matrix J

    
    rLeftBoundary = rLeft(1:end);
% 
%     rGrid = reshape(r, [N_phi, N_theta]);
%     r_eq_temp = 2*rGrid(2,:)-rGrid(1,:);
%     r_mer_temp = 2*rGrid(end,:)-rGrid(end-1,:);
%     
%     r_eq = horzcat(0, r_eq_temp, 0).';
%     r_mer = horzcat(0, r_mer_temp, 0).';
%     
iter=1;
for i=1:N_phi*N_theta
    
    for j=1:N_phi*N_theta
        
        if i==j            
           I(iter)=i;
           J(iter)=j;
           V(iter)=dfkdrk_Numerical_Bottom(i, Max_deg_equ, Max_deg_phi_direction, delta_theta, delta_phi, r, rLeftBoundary, r_eq, r_mer, ThetaLeftBoundaryDeg, SystemParameters, epsilon);
           iter=iter+1;
        
        elseif j==i+1       % sup-diagonal
            I(iter)=i;
            J(iter)=j;            
            V(iter)=dfkdrkp1_Numerical_Bottom(i, Max_deg_equ, Max_deg_phi_direction, delta_theta, delta_phi, r, rLeftBoundary, r_eq, r_mer, ThetaLeftBoundaryDeg, SystemParameters, epsilon);
            iter=iter+1;
            
        elseif j==i-1  % sub-diagonal
            I(iter)=i;
            J(iter)=j;            
            V(iter)= dfkdrkm1_Numerical_Bottom(i, Max_deg_equ, Max_deg_phi_direction, delta_theta, delta_phi, r, rLeftBoundary, r_eq, r_mer, ThetaLeftBoundaryDeg, SystemParameters, epsilon);
            iter=iter+1;            
        
        elseif j==i-N_phi  % N-sub-diagonal
            I(iter)=i;
            J(iter)=j;            
            V(iter)= dfkdrkmN_Numerical_Bottom(i, Max_deg_equ, Max_deg_phi_direction, delta_theta, delta_phi, r, rLeftBoundary, r_eq, r_mer, ThetaLeftBoundaryDeg, SystemParameters, epsilon);
            iter=iter+1;   
     

        elseif j==i+N_phi      % N-sup-diagonal
           I(iter)=i;
           J(iter)=j;                       
           V(iter)=dfkdrkpN_Numerical_Bottom(i, Max_deg_equ, Max_deg_phi_direction, delta_theta, delta_phi, r, rLeftBoundary, r_eq, r_mer, ThetaLeftBoundaryDeg, SystemParameters, epsilon);
           iter=iter+1;
           
        end
        
    end    
end

Jacobian_Numerical_Bottom = sparse(I,J,V);


end