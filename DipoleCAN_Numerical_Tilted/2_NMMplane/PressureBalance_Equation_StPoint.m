function [PressureBalance] = PressureBalance_Equation_StPoint(CoordinateDetails, theta, r, drdtheta, drdphi, ParaSystem)

    Beta = ParaSystem.Beta;

%--- Scaling r to r0
    if strcmp(CoordinateDetails.Units, 'Rp')
        ConvFactor = ParaSystem.Rp ./ ParaSystem.r0;
    elseif strcmp(CoordinateDetails.Units, 'r0')
        ConvFactor = 1;
    end

    r = r*ConvFactor;

% %--- From StPoint-centered coordinates (theta, r) to KSM-framed coordinates (X, Y, Z)
    StPoint_r0 = ParaSystem.Nose.KSMU_Rp * ParaSystem.Rp / ParaSystem.r0;
    StPointDirection = StPoint_r0 ./ norm(StPoint_r0);
    ex = [1; 0; 0];
    AlphaStPoint = acos(dot(ex, StPointDirection));

    %% TEST: Position in StPoint frame
%     Theta_XKSM = theta + AlphaStPoint;
%     XKSM = r*cos(Theta_XKSM);
%     ZKSM = -r*sin(Theta_XKSM);
%     YKSM = 0;
    
%     [~, ZKSM, XKSM] = sph2cart(0, pi/2-theta, r);
%     XKSM = r*cos(theta);
%     ZKSM = r*sin(theta);
%     YKSM = 0;
    
    %% TEST
    
    
    %% Conversion: KSM-U -> StPoint Frame
    Nose = ParaSystem.Nose;
    SSN_Direction = Nose.KSMU_Rp ./ norm (Nose.KSMU_Rp);
    ex_SSN = SSN_Direction;

    MAligned = ParaSystem.M .* [0; 0; 1]; 
    MTilted = ParaSystem.Tilt.RotationPhi*( ParaSystem.Tilt.RotationTheta*MAligned);
    ey_SSN = cross(MTilted, ex_SSN) ./ norm(cross(MTilted, ex_SSN));

    ez_SSN = cross(ex_SSN, ey_SSN);
    M_StPoint2KSMU = [ex_SSN, ey_SSN, ez_SSN];
    
    M_KSMU2StPoint =  (M_StPoint2KSMU)^-1;
    
    %% Position Vector, in KSMU Frame
    X = r*cos(theta);
    Z = r*sin(theta);
    PositionVector_StPoint = [X; 0; Z];
    PositionVector_KSMU = M_StPoint2KSMU * PositionVector_StPoint;
    

%--- Computing Bfield using KSM coordinates
    CoordinateDetails.Type = 'Cartesian_KSM';
    CoordinateDetails.Units = 'r0';
    B_Cartesian_KSMU = BField(CoordinateDetails, PositionVector_KSMU, ParaSystem);   
    B_Cartesian_StPoint = M_KSMU2StPoint*B_Cartesian_KSMU;
    
    CoordinateDetails.Type = 'Spherical_MB';

%--- Solar Wind Velocity Vector, KSM    
    vSW_Cartesian_KSMU = [ -1; 0; 0];
    vSW_Cartesian_StPoint = M_KSMU2StPoint*vSW_Cartesian_KSMU;
    
    
%--- Normal Vector in StPoint-centered Cartesian Coordinates
%     if theta == 0

        % In StPoint directions   
%         n_Cartesian = -vSW_Cartesian_StPoint.';
%         Bfield_Cartesian = B_Cartesian_StPoint;  
%         vSW_Cartesian = vSW_Cartesian_StPoint;

%     else
        Theta = theta;
        Phi = pi/2;
        
        P_cart2sph = [ sin(Theta)*cos(Phi)     sin(Theta)*sin(Phi)    cos(Theta)    ;...
                        cos(Theta)*cos(Phi)     cos(Theta)*sin(Phi)    -sin(Theta)   ;...
                        -sin(Phi)               cos(Phi)               0             ...
                    ];
        n_vec = [ 1 ; (-1/r)*drdtheta ; -1/(r*(sin(Theta)))*drdphi ].';
        if theta == 0
            n_vec(3) = 0;
        end
        
        n_Cartesian_MB = (P_cart2sph.') * n_vec.' ;
        n_Cartesian_StPoint = [n_Cartesian_MB(3), n_Cartesian_MB(1), n_Cartesian_MB(2)];
        n_Cartesian_StPoint =  n_Cartesian_StPoint ./ norm(n_Cartesian_StPoint);

        % In SSN-centered directions   
        n_Cartesian = n_Cartesian_StPoint;
        Bfield_Cartesian = B_Cartesian_StPoint;
        vSW_Cartesian = vSW_Cartesian_StPoint;

%         n_Cartesian_KSMU = M_StPoint2KSMU*n_Cartesian_StPoint.'
%         B_Cartesian_KSMU = M_StPoint2KSMU*B_Cartesian_StPoint

%     end


%--- Setting Up Pressure Balance
    PressureMag =  norm(cross(n_Cartesian, Bfield_Cartesian))*sqrt(1+Beta);
    PressureSW = -(1/2)*dot(n_Cartesian, vSW_Cartesian);
    PressureDifference = PressureMag - PressureSW;

    PressureAverage = (1/2)*(PressureMag+PressureSW);

    % PressureBalance = PressureDifference / (PressureAverage);
    PressureBalance = PressureDifference;
    

end