function [dfkdrk, dfkdrkp1, dfkdrkm1, dfkdrkpN, dfkdrkmN] = Derivatives(k, r, rss, rEqu, rMer, DeltaThetaDeg, DeltaPhiDeg, ThetaMaxDeg, PhiMaxDeg)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Computes the relevant partial derivatives relative to point k %%%%%%
%
% Inputs
%    - DeltaThetaDeg, DeltaPhiDeg: scalars, angular increments, in degrees
%    - ThetaCusp: scalar, theta value of the cusp, in degrees
%    - ThetaMaxDeg, PhiMaxDeg: scalars, maximum values for theta and phi,
%           in degrees
%    - rEquator: vector, solution in the equatorial plane
%           theta = 0:DeltaThetaDeg:ThetaMaxDeg
%    - rMeridian: vector, solution in the noon-midnight meridian plane
%           theta = 0:DeltaThetaDeg:ThetaMaxDeg
%    - SurfaceTot: phi*theta array, surface to be corrected
%           on the entire grid
%    - NumIterationsTop, NumIterationsBottom: scalars, number of iterations
%           for the correction of the top/bottom sub-grid
%
% Outputs
%     - SurfaceCorrected: phi*theta array, corrected surface
%           on the entire grid
%     - rCroppedCorrected: vector, corrected upper sub-grid
%     - SurfaceCroppedExtrapolated: phi*theta array, corrected surface on
%     the sub-upper grid
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% dfkdrk = zeros(length(r), 1);
% dfkdrkp1 = zeros(length(r), 1);
% dfkdrkm1 = zeros(length(r), 1);
% dfkdrkpN = zeros(length(r), 1);
% dfkdrkmN = zeros(length(r), 1);

    Bp = 20000*10^(-9);     % Equatorial field, T
    Rp = 60280*10^3;        % Planet radius, m
    
    mu_0 = 4*pi*10^(-7);
    Pswnpa = 0.02;
    
    b0 = sqrt(2*mu_0*Pswnpa*10^(-9));   % Field scale
    r0 = (2*Bp*Rp^3/b0)^(1/3);          % Distance scale
    M = (Bp/b0)*(Rp/r0)^3;

    Nb_points_theta = ThetaMaxDeg/DeltaThetaDeg+1;    % Number of points on each direction of the grid, from 0 to Max_deg_equ
    N_theta = Nb_points_theta-2; % Number of points on theta-direction of the grid, from 1 to Max_deg_equ-1

    Nb_points_phi= PhiMaxDeg/DeltaPhiDeg+1;    % Number of points on each direction of the grid, from 0 to Max_deg_equ
    N_phi = Nb_points_phi-2; % Number of points on theta-direction of the grid, from 1 to Max_deg_equ-1

    
    
    Max_deg_theta = ThetaMaxDeg-DeltaThetaDeg;
    Max_deg_phi= PhiMaxDeg-DeltaPhiDeg;
    
    delta_theta = DeltaThetaDeg*(pi/180);
    delta_phi = DeltaPhiDeg*(pi/180);
    
    theta_vec_1 = repmat(DeltaThetaDeg:DeltaThetaDeg:Max_deg_theta,[N_phi 1]);
    theta_vec = theta_vec_1(:)*pi/180;

    phi_vec_1 = repmat(DeltaPhiDeg:DeltaPhiDeg:Max_deg_phi,[1 N_theta]).';
    phi_vec = phi_vec_1(:)*pi/180;
    


syms rk rkp1 rkm1 rkpN rkmN;

    if k == 1
        
        n_vec = [ 1, -1/rk*(rkpN-rss) / (2*delta_theta), -1/(rk*sin(theta_vec(k))) * (rkp1-rEqu(2)) / (2*delta_phi)];
        B = M*(1/rk)^3 * [ -2*sin(theta_vec(k))*sin(phi_vec(k)), cos(theta_vec(k))*sin(phi_vec(k)), cos(phi_vec(k)) ];
        v = [ -cos(theta_vec(k)), sin(theta_vec(k)), 0 ];
        PressureBalance =   norm(cross(n_vec, B)) + (1/2)*dot(n_vec,v);  
        
        dfkdrkFunction =   matlabFunction(diff(PressureBalance, rk));
        dfkdrkp1Function = matlabFunction(diff(PressureBalance, rkp1));
        dfkdrkpNFunction = matlabFunction(diff(PressureBalance, rkpN));
        
        dfkdrk = dfkdrkFunction(r(k), r(k+1), r(k+N_phi));
        dfkdrkp1 = dfkdrkp1Function(r(k), r(k+1), r(k+N_phi));
        dfkdrkm1 = 0;
        dfkdrkpN = dfkdrkpNFunction(r(k), r(k+1), r(k+N_phi));
        dfkdrkmN = 0;
        
        
    elseif k == N_phi
        
        n_vec = [ 1, -1/rk*(rkpN-rss) / (2*delta_theta), -1/(rk*sin(theta_vec(k))) * (rMer(2)-rkm1) / (2*delta_phi)];
        B = M*(1/rk)^3 * [ -2*sin(theta_vec(k))*sin(phi_vec(k)), cos(theta_vec(k))*sin(phi_vec(k)), cos(phi_vec(k)) ];
        v = [ -cos(theta_vec(k)), sin(theta_vec(k)), 0 ];
        PressureBalance =   norm(cross(n_vec, B)) + (1/2)*dot(n_vec,v);  
        
        dfkdrkFunction =   matlabFunction(diff(PressureBalance, rk));
        dfkdrkm1Function = matlabFunction(diff(PressureBalance, rkm1));
        dfkdrkpNFunction = matlabFunction(diff(PressureBalance, rkpN));
                
        dfkdrk = dfkdrkFunction(r(k), r(k-1), r(k+N_phi) );
        dfkdrkp1 = 0;
        dfkdrkm1 = dfkdrkm1Function(r(k), r(k-1), r(k+N_phi));
        dfkdrkpN = dfkdrkpNFunction(r(k), r(k-1), r(k+N_phi));
        dfkdrkmN = 0;

    elseif k == N_phi*N_theta-N_phi+1

        n_vec = [ 1, -1/rk*(2*rk-rkmN-rkmN)/ (2*delta_theta), -1/(rk*sin(theta_vec(k))) * (rkp1-rEqu(N_theta+1)) / (2*delta_phi)];
        B = M*(1/rk)^3 * [ -2*sin(theta_vec(k))*sin(phi_vec(k)), cos(theta_vec(k))*sin(phi_vec(k)), cos(phi_vec(k)) ];
        v = [ -cos(theta_vec(k)), sin(theta_vec(k)), 0 ];    
        PressureBalance =   norm(cross(n_vec, B)) + (1/2)*dot(n_vec,v);  
        
        dfkdrkFunction =   matlabFunction(diff(PressureBalance, rk));
        dfkdrkp1Function = matlabFunction(diff(PressureBalance, rkp1));
        dfkdrkmNFunction = matlabFunction(diff(PressureBalance, rkmN));
        
        dfkdrk = dfkdrkFunction(r(k), r(k+1), r(k-N_phi));
        dfkdrkp1 = dfkdrkp1Function(r(k), r(k+1), r(k-N_phi));
        dfkdrkm1 = 0;
        dfkdrkpN = 0;
        dfkdrkmN = dfkdrkmNFunction(r(k), r(k+1), r(k-N_phi));
        
    elseif k == N_phi*N_theta
        
        n_vec = [ 1, -1/rk*(2*rk-rkmN-rkmN)/ (2*delta_theta), -1/(rk*sin(theta_vec(k))) * (rMer(N_theta+1)-rkm1) / (2*delta_phi)];
        B = M*(1/rk)^3 * [ -2*sin(theta_vec(k))*sin(phi_vec(k)), cos(theta_vec(k))*sin(phi_vec(k)), cos(phi_vec(k)) ];
        v = [ -cos(theta_vec(k)), sin(theta_vec(k)), 0 ];      
        PressureBalance =   norm(cross(n_vec, B)) + (1/2)*dot(n_vec,v);  
        
        dfkdrkFunction = matlabFunction(diff(PressureBalance, rk));
        dfkdrkm1Function = matlabFunction(diff(PressureBalance, rkm1));
        dfkdrkmNFunction = matlabFunction(diff(PressureBalance, rkmN));
        
        dfkdrk = dfkdrkFunction(r(k), r(k-1), r(k-N_phi));
        dfkdrkp1 = 0;
        dfkdrkm1 = dfkdrkm1Function(r(k), r(k-1), r(k-N_phi));
        dfkdrkpN = 0;
        dfkdrkmN = dfkdrkmNFunction(r(k), r(k-1), r(k-N_phi));        
        
    elseif k>1 && k<N_phi   %%% left boundary

        n_vec = [ 1, -1/rk*(rkpN-rss) / (2*delta_theta), -1/(rk*sin(theta_vec(k))) * (rkp1-rkm1) / (2*delta_phi)];
        B = M*(1/rk)^3 * [ -2*sin(theta_vec(k))*sin(phi_vec(k)), cos(theta_vec(k))*sin(phi_vec(k)), cos(phi_vec(k)) ];
        v = [ -cos(theta_vec(k)), sin(theta_vec(k)), 0 ];
        PressureBalance =   norm(cross(n_vec, B)) + (1/2)*dot(n_vec,v);  
        
        dfkdrkFunction = matlabFunction(diff(PressureBalance, rk));
        dfkdrkp1Function = matlabFunction(diff(PressureBalance, rkp1));
        dfkdrkm1Function = matlabFunction(diff(PressureBalance, rkm1));
        dfkdrkpNFunction = matlabFunction(diff(PressureBalance, rkpN));
        
        dfkdrk = dfkdrkFunction(r(k), r(k+1), r(k-1), r(k+N_phi));
        dfkdrkp1 = dfkdrkp1Function(r(k), r(k+1), r(k-1), r(k+N_phi));
        dfkdrkm1 = dfkdrkm1Function(r(k), r(k+1), r(k-1), r(k+N_phi));
        dfkdrkpN = dfkdrkpNFunction(r(k), r(k+1), r(k-1), r(k+N_phi));
        dfkdrkmN = 0;
        
    elseif k>1 && mod(k,N_phi)==1 && k<N_phi*N_theta-N_phi+1    %%% bottom boundary

        n_vec = [ 1, -1/rk*(rkpN-rkmN) / (2*delta_theta), -1/(rk*sin(theta_vec(k))) * (rkp1-rEqu((k-1)/N_phi+2)) / (2*delta_phi)];
        B = M*(1/rk)^3 * [ -2*sin(theta_vec(k))*sin(phi_vec(k)), cos(theta_vec(k))*sin(phi_vec(k)), cos(phi_vec(k)) ];
        v = [ -cos(theta_vec(k)), sin(theta_vec(k)), 0 ];
        PressureBalance =   norm(cross(n_vec, B)) + (1/2)*dot(n_vec,v);  
        
        dfkdrkFunction = matlabFunction(diff(PressureBalance, rk));
        dfkdrkp1Function = matlabFunction(diff(PressureBalance, rkp1));
        dfkdrkpNFunction = matlabFunction(diff(PressureBalance, rkpN));
        dfkdrkmNFunction = matlabFunction(diff(PressureBalance, rkmN));
        
        dfkdrk = dfkdrkFunction(r(k), r(k+1), r(k+N_phi), r(k-N_phi));
        dfkdrkp1 = dfkdrkp1Function(r(k), r(k+1), r(k+N_phi), r(k-N_phi));
        dfkdrkm1 = 0;
        dfkdrkpN = dfkdrkpNFunction(r(k), r(k+1), r(k+N_phi), r(k-N_phi));
        dfkdrkmN = dfkdrkmNFunction(r(k), r(k+1), r(k+N_phi), r(k-N_phi));        
                
    elseif k>N_phi && mod(k,N_phi)==0 && k<N_phi*N_theta  %%% upper boundary
        
        n_vec = [ 1, -1/rk*(rkpN-rkmN) / (2*delta_theta), -1/(rk*sin(theta_vec(k))) * (rMer(k/N_phi+1)-rkm1) / (2*delta_phi)];
        B = M*(1/rk)^3 * [ -2*sin(theta_vec(k))*sin(phi_vec(k)), cos(theta_vec(k))*sin(phi_vec(k)), cos(phi_vec(k)) ];
        v = [ -cos(theta_vec(k)), sin(theta_vec(k)), 0 ];    
        PressureBalance =   norm(cross(n_vec, B)) + (1/2)*dot(n_vec,v);  
        
        dfkdrkFunction = matlabFunction(diff(PressureBalance, rk));
        dfkdrkm1Function = matlabFunction(diff(PressureBalance, rkm1));
        dfkdrkpNFunction = matlabFunction(diff(PressureBalance, rkpN));
        dfkdrkmNFunction = matlabFunction(diff(PressureBalance, rkmN));
        
        dfkdrk = dfkdrkFunction(r(k), r(k-1), r(k+N_phi), r(k-N_phi));
        dfkdrkp1 = 0;
        dfkdrkm1 = dfkdrkm1Function(r(k), r(k-1), r(k+N_phi), r(k-N_phi));
        dfkdrkpN = dfkdrkpNFunction(r(k), r(k-1), r(k+N_phi), r(k-N_phi));
        dfkdrkmN = dfkdrkmNFunction(r(k), r(k-1), r(k+N_phi), r(k-N_phi));        
                        
    elseif k>N_phi*N_theta-N_phi+1 && k<N_phi*N_theta  %%% right boundary

        n_vec = [ 1, -1/rk*(2*rk-rkmN-rkmN) / (2*delta_theta), -1/(rk*sin(theta_vec(k))) * (rkp1-rkm1) / (2*delta_phi)];
        B = M*(1/rk)^3 * [ -2*sin(theta_vec(k))*sin(phi_vec(k)), cos(theta_vec(k))*sin(phi_vec(k)), cos(phi_vec(k)) ];
        v = [ -cos(theta_vec(k)), sin(theta_vec(k)), 0 ];
        PressureBalance =   norm(cross(n_vec, B)) + (1/2)*dot(n_vec,v);  
        
        dfkdrkFunction = matlabFunction(diff(PressureBalance, rk));
        dfkdrkp1Function = matlabFunction(diff(PressureBalance, rkp1));
        dfkdrkm1Function = matlabFunction(diff(PressureBalance, rkm1));
        dfkdrkmNFunction = matlabFunction(diff(PressureBalance, rkmN));
        
        dfkdrk = dfkdrkFunction(r(k), r(k+1), r(k-1), r(k-N_phi));
        dfkdrkp1 = dfkdrkp1Function(r(k), r(k+1), r(k-1), r(k-N_phi));
        dfkdrkm1 = dfkdrkm1Function(r(k), r(k+1), r(k-1), r(k-N_phi));
        dfkdrkpN = 0;
        dfkdrkmN = dfkdrkmNFunction(r(k), r(k+1), r(k-1), r(k-N_phi));              
        
    else

        n_vec = [ 1, -1/rk*(rkpN-rkmN) / (2*delta_theta), -1/(rk*sin(theta_vec(k))) * (rkp1-rkm1) / (2*delta_phi)];
        B = M*(1/rk)^3 * [ -2*sin(theta_vec(k))*sin(phi_vec(k)), cos(theta_vec(k))*sin(phi_vec(k)), cos(phi_vec(k)) ];
        v = [ -cos(theta_vec(k)), sin(theta_vec(k)), 0 ];
        PressureBalance =   norm(cross(n_vec, B)) + (1/2)*dot(n_vec,v);  
        
        dfkdrkFunction = matlabFunction(diff(PressureBalance, rk));
        dfkdrkp1Function = matlabFunction(diff(PressureBalance, rkp1));
        dfkdrkm1Function = matlabFunction(diff(PressureBalance, rkm1));
        dfkdrkpNFunction = matlabFunction(diff(PressureBalance, rkpN));
        dfkdrkmNFunction = matlabFunction(diff(PressureBalance, rkmN));
        
        dfkdrk = dfkdrkFunction(r(k), r(k+1), r(k-1), r(k+N_phi), r(k-N_phi));
        dfkdrkp1 = dfkdrkp1Function(r(k), r(k+1), r(k-1), r(k+N_phi), r(k-N_phi));
        dfkdrkm1 = dfkdrkm1Function(r(k), r(k+1), r(k-1), r(k+N_phi), r(k-N_phi));
        dfkdrkpN = dfkdrkpNFunction(r(k), r(k+1), r(k-1), r(k+N_phi), r(k-N_phi));
        dfkdrkmN = dfkdrkmNFunction(r(k), r(k+1), r(k-1), r(k+N_phi), r(k-N_phi));
        
    end
    
end