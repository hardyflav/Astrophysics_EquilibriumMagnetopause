function PressureBalanceFringingFieldv2 = PressureBalanceFringingFieldv2( k, r, rss, rEquator, rRightBoundary, rMeridian, DeltaThetaDeg, DeltaPhiDeg, ThetaMaxDeg, PhiMaxDeg, BFringingThetaInterp, BFringingPhiInterp, BFringingrInterp, BIMF )
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs
%     - k: scalar, position on the grid described by a linear array
%     - r: phi*theta INTERIOR grid, surface described by a linear array
%     - rss: scalar, value of the stand-off distance for the left boundary
%     - rEquator: vector, equatorial solution for the bottom boundary
%           length: ThetaMaxDeg/DeltaThetaDeg+1
%           theta = 0 : DeltaThetaDeg : MaxThetaDeg
%     - rMeridian: vector, meridional solution for the top boundary
%           length: ThetaMaxDeg/DeltaThetaDeg+1
%           theta = 0 : DeltaThetaDeg : MaxThetaDeg
%     - rRightBoundary: vector, linear expansion for the right boundary
%           length: ThetaMaxDeg/DeltaThetaDeg+1
%           theta = 0 : DeltaThetaDeg : MaxThetaDeg
%     - DeltaThetaDeg, DeltaPhiDeg: Scalars, angular increments in degrees
%     - BFringingThetaInterp: phi*theta INTERIOR grid of Bf along etheta
%     - BFringingPhiInterp: phi*theta INTERIOR grid of Bf along ephi
%     - BFringingrInterp: phi*theta INTERIOR grid of Bf along er
%     - BIMF: 1x3 vector, value of the external interplanetary field,
%           in CARTESIAN coordinate system [Bx, By, Bz] 
%
%
% Output
%     - PressureBalanceFringingFieldv2: scalar, pressure balance at point k
%    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%% Planetary constants

    Bp = 20000*10^(-9);     % Equatorial field, T
    Rp = 60280*10^3;        % Planet radius, m
    
    mu_0 = 4*pi*10^(-7);
    Pswnpa = 0.02;
    
    b0 = sqrt(2*mu_0*Pswnpa*10^(-9));   % Field scale
    r0 = (2*Bp*Rp^3/b0)^(1/3);          % Distance scale
    M = (Bp/b0)*(Rp/r0)^3;

%%% Definition of the [theta*phi] grid

    Nb_points_theta = ThetaMaxDeg/DeltaThetaDeg+1;    % Number of points on each direction of the grid, from 0 to ThetaMaxDeg
    N_theta = Nb_points_theta-2; % Number of points on theta-direction of the grid, from 1 to ThetaMaxDeg-1
    Nb_points_phi= PhiMaxDeg/DeltaPhiDeg+1;    % Number of points on each direction of the grid, from 0 to ThetaMaxDeg
    N_phi = Nb_points_phi-2; % Number of points on theta-direction of the grid, from 1 to ThetaMaxDeg-1

    Max_deg_theta = ThetaMaxDeg-DeltaThetaDeg;
    Max_deg_phi= PhiMaxDeg-DeltaPhiDeg;
    delta_theta = DeltaThetaDeg*(pi/180);
    delta_phi = DeltaPhiDeg*(pi/180);
    
    theta_vec_1 = repmat(DeltaThetaDeg:DeltaThetaDeg:Max_deg_theta,[N_phi 1]);
    theta_vec = theta_vec_1(:)*pi/180;

    phi_vec_1 = repmat(DeltaPhiDeg:DeltaPhiDeg:Max_deg_phi,[1 N_theta]).';
    phi_vec = phi_vec_1(:)*pi/180;


%%% Boundary conditions on the interior grid

    rLeft = rss;
    rBottom = rEquator(2:end-1);
    rTop = rMeridian(2:end-1);
    rRight = rRightBoundary(2:end-1);

%%% Extraction of the Components for the fringing field BFringing, as linear
%%% arrays
    
    BfTheta = BFringingThetaInterp(:);
    BfPhi = BFringingPhiInterp(:);
    Bfr = BFringingrInterp(:);
    
    Bfk = [Bfr(k), BfTheta(k), BfPhi(k)];
    
    MatrixCartesianToSpherical = [sin(theta_vec(k))*cos(phi_vec(k)) cos(theta_vec(k))*cos(phi_vec(k)) -sin(phi_vec(k)); sin(theta_vec(k))*sin(phi_vec(k)) cos(theta_vec(k))*sin(phi_vec(k)) cos(phi_vec(k)); cos(theta_vec(k)) -sin(theta_vec(k)) 0];
    BIMFSpherical = BIMF*MatrixCartesianToSpherical;
    
%%% Back to the pressure balance at r(k) = P, on the surface    
    
        if k == 1
            n_vec = [ 1, -1/r(k)*(r(k+N_phi)-rLeft) / (2*delta_theta), -1/(r(k)*sin(theta_vec(k))) * (r(k+1)-rBottom(1)) / (2*delta_phi)];
        
        elseif k == N_phi
            n_vec = [ 1, -1/r(k)*(r(k+N_phi)-rLeft) / (2*delta_theta), -1/(r(k)*sin(theta_vec(k))) * (rTop(1)-r(k-1)) / (2*delta_phi)];
        
        elseif k == N_phi*N_theta-N_phi+1
            n_vec = [ 1, -1/r(k)*(rRight(1) -r(k-N_phi))/ (2*delta_theta), -1/(r(k)*sin(theta_vec(k))) * (r(k+1)-rBottom(end)) / (2*delta_phi)];
            
        elseif k == N_phi*N_theta
            n_vec = [ 1, -1/r(k)*(rRight(end)-r(k-N_phi))/ (2*delta_theta), -1/(r(k)*sin(theta_vec(k))) * (rTop(end)-r(k-1)) / (2*delta_phi)];
            
        elseif k>1 && k<N_phi   % left boundary
        	n_vec = [ 1, -1/r(k)*(r(k+N_phi)-rLeft) / (2*delta_theta), -1/(r(k)*sin(theta_vec(k))) * (r(k+1)-r(k-1)) / (2*delta_phi)];
            
        elseif k>1 && mod(k,N_phi)==1 && k<N_phi*N_theta-N_phi+1    % bottom boundary
            n_vec = [ 1, -1/r(k)*(r(k+N_phi)-r(k-N_phi)) / (2*delta_theta), -1/(r(k)*sin(theta_vec(k))) * (r(k+1)-rBottom((k-1)/N_phi+1)) / (2*delta_phi)];
            
        elseif k>N_phi && mod(k,N_phi)==0 && k<N_phi*N_theta  % upper boundary
            n_vec = [ 1, -1/r(k)*(r(k+N_phi)-r(k-N_phi)) / (2*delta_theta), -1/(r(k)*sin(theta_vec(k))) * (rTop(k/N_phi)-r(k-1)) / (2*delta_phi)];
                
        elseif k>N_phi*N_theta-N_phi+1 && k<N_phi*N_theta  % right boundary
            n_vec = [ 1, -1/r(k)*(rRight( k-(N_phi*N_theta-N_phi) )-r(k-N_phi)) / (2*delta_theta), -1/(r(k)*sin(theta_vec(k))) * (r(k+1)-r(k-1)) / (2*delta_phi)];
   
        else
            n_vec = [ 1, -1/r(k)*(r(k+N_phi)-r(k-N_phi)) / (2*delta_theta), -1/(r(k)*sin(theta_vec(k))) * (r(k+1)-r(k-1)) / (2*delta_phi)];
        
        end

     Bk = M*(1/r(k))^3 * [ -2*sin(theta_vec(k))*sin(phi_vec(k)), cos(theta_vec(k))*sin(phi_vec(k)), cos(phi_vec(k)) ];
     Btot = Bk+Bfk.';
     v = [ -cos(theta_vec(k)), sin(theta_vec(k)), 0 ];        

    PressureBalanceFringingFieldv2 = norm(cross(n_vec, Btot - BIMFSpherical)) + (1/2)*dot(n_vec,v);

end