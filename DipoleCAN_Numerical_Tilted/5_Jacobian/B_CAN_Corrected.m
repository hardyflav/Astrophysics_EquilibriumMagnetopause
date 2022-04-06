function B_CAN_Cartesian = B_CAN_Corrected(CAN_DiskParameters, r, theta, phi)
%
% provides cylindrical components of field for axisymmetric ring current
% extending from a to b, with half-thickness D, and current parameter I
%
% input: x = [I a b D]
%        xdata = [z; rho] in 2 column format (phi is irrelevant)
%
% output: Bcyl_disk, contains both components sequentially
% ie
% Bz_disk = Bcyl_disk(1:length(Bcyl_disk)/2);
% Brho_disk= Bcyl_disk(length(Bcyl_disk)/2+1:length(Bcyl_disk));
%
% Giacomo Giampieri, October 2002


% System parameters
    Bp = 20000*10^(-9);     % Equatorial field, T
    Rp = 60280*10^3;        % Planet radius, m
    
    mu_0 = 4*pi*10^(-7);
    Pswnpa = 0.02;
    
    b0 = sqrt(2*mu_0*Pswnpa*10^(-9));   % Field scale
    r0 = (2*Bp*Rp^3/b0)^(1/3);          % Distance scale
    M = (Bp/b0)*(Rp/r0)^3;
    
    z = 0.* r;
    rho = 0.* r;
    
% Converting the position from MB64 to CAN89, in Rp
%     Elevation = wrapToPi( theta + pi/2 );
    Elevation = ( pi/2 - theta );
   
    [Y_KSM, Z_KSM, X_KSM] = sph2cart(phi, Elevation, r);
%     [Y_KSM, Z_KSM, X_KSM] = sph2cart(phi, pi/2-theta, r);
    z = Z_KSM;
    rho = sqrt( X_KSM.^2 + Y_KSM.^2 );
    
%     Cond1 = (theta == 0);
% 	rho(Cond1) = r(Cond1) ;
% 
%     Cond2 = (theta == pi/2 & phi == pi/2);
%     z(Cond2) = r(Cond2) ;
%     
%     Cond3 = (theta == pi/2 & phi == -pi/2);
%     z(Cond3) = -r(Cond3) ;
    
% Converting the position from MB64 to CAN89, in Rp
%     if theta == 0
%         z = 0;
%         rho = r *r0/Rp;
%     elseif theta == pi/2 && phi == pi/2
%         z = r *r0/Rp;
%         rho = 0;
%     else
%          z = r.*sin(theta).*sin(phi)  .*r0/Rp ;
%         rho = r.*sqrt( 1-(sin(theta).*sin(phi)).^2 )  .*r0/Rp;
%     end
    
% Computing the CAN disk field components in CAN89
    x = CAN_DiskParameters;
    xdata = [z(:).'; rho(:).'];
    B_CAN = diskfield_cyl(x, xdata);
    
    B_CAN_Z = B_CAN(1:length(B_CAN)/2) .';    
    B_CAN_rho = B_CAN(length(B_CAN)/2+1:length(B_CAN)) .';
    
    [rowNaN] = find( isnan(B_CAN_Z) ) ;
    B_CAN_Z(rowNaN) = (1/2) .* ( B_CAN_Z(rowNaN-1) + B_CAN_Z(rowNaN+1) );
    
    B_CAN_rho(theta==pi/2) = 0;
            
    % Adjusting fields at singular points
%     rhoVec = rho(:);
%     phiVec = phi(:);
%     thetaVec = theta(:);
%     B_CAN_rho(rhoVec == 0) = 0;
%     B_CAN_rho(phiVec == 0) = 0;
%     
%     B_CAN_Z( isnan(B_CAN_Z) ) = 0;
    
    %{
    rhoVec( phiVec==0 & thetaVec == pi/2)
    r( phiVec==0 & thetaVec == pi/2)
    
    phiVec( isnan(B_CAN_rho))
    thetaVec( isnan(B_CAN_rho))
    rVec = r(:);
    B_CAN_rho( isnan(B_CAN_rho))
    zVec(thetaVec == pi/2)
    
    B_CAN_Z(phiVec == 0) = 0;
    %}
    
    B_CAN_ZGrid = reshape(B_CAN_Z, size(z(:)));
    B_CAN_rhoGrid = reshape(B_CAN_rho, size(z(:)));
        
% Converting the CAN disk field components in MB64 Coordinate system
    
    B_CAN_r = 0.* B_CAN_ZGrid;
    B_CAN_theta = 0.* B_CAN_ZGrid;
    B_CAN_phi = 0.* B_CAN_ZGrid;
    
    [XMB, YMB, ZMB] = sph2cart(phi, pi/2-theta, r);
    XKSM = ZMB;
    YKSM = XMB;
    ZKSM = YMB;
    
    [Phi_XY_KSM, ~, ~] = cart2sph(XKSM, YKSM, ZKSM);
    
    
%     PhiList = linspace(0, 2*pi, length(z(:))) .';
    PhiVec = Phi_XY_KSM(:);
    B_CAN_X = cos(PhiVec) .* B_CAN_rhoGrid - sin(PhiVec) .* B_CAN_phi ;
    B_CAN_Y = sin(PhiVec) .* B_CAN_rhoGrid + cos(PhiVec) .* B_CAN_phi ;
    
%     B_CAN_X( abs(phi) == pi/2 ) = 0;
%     B_CAN_Y( abs(phi) == pi/2 ) = 0;
%     B_CAN_Z( (phi) == -pi/2 ) = - B_CAN_Z( (phi) == pi/2 );
    
    B_CAN_Cartesian.X = reshape(B_CAN_X, size(r));
    B_CAN_Cartesian.Y = reshape(B_CAN_Y, size(r));
    B_CAN_Cartesian.Z = reshape(B_CAN_ZGrid, size(r));
    
    
    %{
	[X_CAN, Y_CAN, Z_CAN] = pol2cart(0.*B_CAN_rhoGrid, B_CAN_rhoGrid, B_CAN_ZGrid);
    
	X_MB = Y_CAN;
	Y_MB = Z_CAN;
	Z_MB = X_CAN;
    
	[B_CAN_phi, ~, B_CAN_r] = cart2sph(X_MB, Y_MB, Z_MB);

        B_CAN_theta = atan2(sqrt(X_MB^2+ Y_MB^2), Z_MB);
    
    Cond3 = (theta == 0 & phi == 0);
    B_CAN_theta(Cond3) = B_CAN_ZGrid(Cond3);
    
    Cond4 = (theta == pi/2 & phi == pi/2);
    B_CAN_r(Cond4) = B_CAN_ZGrid(Cond4);
    
    Cond5 = (phi == 0);
    B_CAN_phi(Cond5) = B_CAN_ZGrid(Cond5);   
    
    Cond6 = (phi == pi/2);
    B_CAN_r(Cond6) = B_CAN_rhoGrid(Cond6) .* cos(theta(Cond6)) + B_CAN_ZGrid(Cond6) .* sin(theta(Cond6));
    B_CAN_theta(Cond6) = -B_CAN_rhoGrid(Cond6) .* sin(theta(Cond6)) + B_CAN_ZGrid(Cond6) .* cos(theta(Cond6));

        %     B_CAN_z = B_CAN(1);
%     B_CAN_rho = B_CAN(2);

% Converting the CAN disk field components in MB64 Coordinate system
    if theta == 0 && phi == 0
        B_CAN_r = 0;
        B_CAN_theta = B_CAN_z;
        B_CAN_phi = 0;   
    elseif theta == pi/2 && phi == pi/2
        B_CAN_r = B_CAN_z;
        B_CAN_theta = 0;
        B_CAN_phi = 0;  
    elseif phi == 0
        B_CAN_r = 0;
        B_CAN_theta = 0;
        B_CAN_phi = B_CAN_z;   
    elseif phi == pi/2
        B_CAN_r = B_CAN_rho*cos(theta) + B_CAN_z*sin(theta);
        B_CAN_theta = -B_CAN_rho*sin(theta) + B_CAN_z*cos(theta);
        B_CAN_phi = 0;   
% 
%     elseif theta == pi/2 && phi == pi/2
%         B_CAN_r = B_CAN_z;
%         B_CAN_theta = 0;
%         B_CAN_phi = 0;   
% %     elseif theta <= 0*pi/180
% %         B_CAN_r = 0;
% %         B_CAN_theta = 0;
%         B_CAN_phi = B_CAN_z;
%         
%     elseif phi == 0
%         B_CAN_r = 0;
%         B_CAN_theta = B_CAN_z;
%         B_CAN_phi = 0;  
%         
%     elseif phi == pi/2
%         X_MB = 0;
%         Y_MB = B_CAN_z;
%         Z_MB = B_CAN_rho;
%         
% %         [B_CAN_phi, ~, B_CAN_r] = cart2sph(X_MB, Y_MB, Z_MB);
%         [B_CAN_phi, B_CAN_r] = cart2pol(Z_MB, Y_MB);
%         B_CAN_theta = 0;
% %         B_CAN_theta = atan2(sqrt(X_MB^2+ Y_MB^2), Z_MB);        
    else
        [X_CAN, Y_CAN, Z_CAN] = pol2cart(0, B_CAN_rho, B_CAN_z);
    
        X_MB = Y_CAN;
        Y_MB = Z_CAN;
        Z_MB = X_CAN;
    
        [B_CAN_phi, ~, B_CAN_r] = cart2sph(X_MB, Y_MB, Z_MB);

        B_CAN_theta = atan2(sqrt(X_MB^2+ Y_MB^2), Z_MB);
    end

    B_CAN_Corrected = [B_CAN_r; B_CAN_theta; B_CAN_phi];
%}

% XPlot = z(:);
% YPlot = rho(:);
% ZPlot = B_CAN_ZGrid;
% scatter3(XPlot, YPlot, ZPlot)

end    
