

function Bf = FringingField(Position, Surface, BfPrevious, BIMF, C, ThetaMaxDeg, PhiMaxDeg, DeltaThetaDeg, DeltaPhiDeg, rss, rEquator, rMeridian, rRightBoundary)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs
%     - Position: 3x1 vector, position at which Bf is calculated (r(M))
%           in SPHERICAL coordinate system [theta, phi, r] 
%     - Surface: Array describing the surface on the inner grid
%     - BfPrevious: 3x1 vector, previous value of the fringing field
%           at Position
%           in Cartesian coordinate system [Bx, By, Bz]
%     - BIMF: 3x1 vector, value of the external interplanetary field,
%           in CARTESIAN coordinate system [Bx, By, Bz] 
%     - C: scalar, 'jumping increment', AnglesNew = AnglesPrevious(1:C:end)
%     - ThetaMaxDeg, PhiMaxDeg: Scalars, maximum values of the angles, on
%           the TOTAL grid (in degrees)
%     - DeltaThetaDeg, DeltaPhiDeg: Scalars, angular increments in degrees
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
%
%
% Output
%     - Bf: 3x1 vector, value of the updated fringing field
%           at Position
%           in Cartesian coordinate system [Bx, By, Bz]
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

    SurfaceLight = Surface(1:C:end, 1:C:end);
    r = SurfaceLight(:);

    Nb_points_theta = ThetaMaxDeg/DeltaThetaDeg+1;
    N_theta = Nb_points_theta-2; % Number of points on theta-direction of the grid, from 1 to ThetaMaxDeg-1
    Nb_points_phi= PhiMaxDeg/DeltaPhiDeg+1;    % Number of points on each direction of the grid, from 0 to ThetaMaxDeg
    N_phi = Nb_points_phi-2; % Number of points on theta-direction of the grid, from 1 to PhiMaxDeg-1

    if (mod(N_theta,C))>0
       NThetaLight = (N_theta-mod(N_theta,C))/C+1;
    else
       NThetaLight = (N_theta-mod(N_theta,C))/C;
    end

    if (mod(N_phi,C))>0
        NPhiLight = (N_phi-mod(N_phi,C))/C+1;
    else
        NPhiLight = (N_phi-mod(N_phi,C))/C;
    end
    N_theta = NThetaLight;
    N_phi = NPhiLight;

    theta_vec_1Light = repmat(DeltaThetaDeg:C*DeltaThetaDeg:ThetaMaxDeg-DeltaThetaDeg,[NPhiLight 1]);
    theta_vecLight = theta_vec_1Light(:)*pi/180;
    phi_vec_1Light = repmat(DeltaPhiDeg:C*DeltaPhiDeg:PhiMaxDeg-DeltaPhiDeg,[1 NThetaLight]).';
    phi_vecLight = phi_vec_1Light(:)*pi/180;
    
    theta_vec = theta_vecLight;
    phi_vec = phi_vecLight;
    
    
    Max_deg_theta = max(theta_vec_1Light(:));
    Max_deg_phi= max(phi_vec_1Light(:));
    delta_theta = DeltaThetaDeg*(pi/180);
    delta_phi = DeltaPhiDeg*(pi/180);
    

%%% Boundary conditions on the interior grid

    rLeft = rss;
    rBottom = rEquator(2:end-1);
    rTop = rMeridian(2:end-1);
    rRight = rRightBoundary(2:end-1);

%%% Iterative construction of Bf - summation over every 'patch' of the
%%% surface

    Bf = [0, 0, 0].'; % Initial value of the updated fringing field
    
    for k = 1 : N_theta*N_phi
        
            
    % (er, etheta, ephi) = SpherToCart * (ex, ey, ez)    
    SpherToCartP = [ sin(theta_vec(k))*cos(phi_vec(k))  sin(theta_vec(k))*sin(phi_vec(k))  cos(theta_vec(k)) ; cos(theta_vec(k))*cos(phi_vec(k))  cos(theta_vec(k))*sin(phi_vec(k))  -sin(theta_vec(k)) ; -sin(phi_vec(k))  cos(phi_vec(k))  0];  
    SpherToCartM = [ sin(Position(1))*cos(Position(2))  sin(Position(1))*sin(Position(2))  cos(Position(1)) ; cos(Position(1))*cos(Position(2))  cos(Position(1))*sin(Position(2))  -sin(Position(1)) ; -sin(Position(2))  cos(Position(2))  0];  

    % Position vectors in Cartesian system
    OPVectorCart = ([r(k), 0, 0]*SpherToCartP).';
    OMVectorCart = ([Position(3), 0, 0]*SpherToCartM).';
    PMVectorCart = - OPVectorCart + OMVectorCart;
    
    if norm(PMVectorCart)>0
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
        
        n_vec_CartP = (n_vec*SpherToCartP).';
    
        % Planetary field at r(k) (on the surface)
        B = M*(1/r(k))^3 * [ -2*sin(theta_vec(k))*sin(phi_vec(k)), cos(theta_vec(k))*sin(phi_vec(k)), cos(phi_vec(k)) ];
        BCartP = (B*SpherToCartP).';
        
        % Previous Fringing field at r(k) (on the surface)
        BfPreviousCartP =  (BfPrevious*SpherToCartP).';
    

        % Contribution of the theta*phi patch around r(k) at P to the fringing field
        % at M
        dBfCartM = 1/(2*pi)*cross( cross(-n_vec_CartP, BCartP + BfPreviousCartP - BIMF), PMVectorCart/norm(PMVectorCart)^3 ) * r(k)^2*sin(theta_vec(k))*delta_theta*delta_phi;
        Bf = Bf + dBfCartM;   
    end   
        
    end
    
end  