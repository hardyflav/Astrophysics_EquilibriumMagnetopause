function [PressureBalance] = PressureBalance_Equation(CoordinateDetails, PointCoordinates, drdtheta_ini, drdphi_ini, ParaSystem)

    Beta = ParaSystem.Beta;
    
%--- Scaling r to r0
    if strcmp(CoordinateDetails.Units, 'Rp')
        ConvFactor = ParaSystem.Rp ./ ParaSystem.r0;
    elseif strcmp(CoordinateDetails.Units, 'r0')
        ConvFactor = 1;
    end
    
    PointCoordinates_KSMU = PointCoordinates .* ConvFactor;
    drdtheta = drdtheta_ini * ConvFactor;
    drdphi = drdphi_ini * ConvFactor;
    CoordinateDetails.Units = 'r0';
    
%--- Solar Wind Velocity Vector, KSMU    
    vSW_Cartesian_KSMU = [ -1; 0; 0];
        
%--- Bfield, KSMU    
    B_Cartesian_KSMU = BField(CoordinateDetails, PointCoordinates_KSMU, ParaSystem);  
    
%--- Normal Vector    

    if CoordinateDetails.SSN == 1
        
    % - In KSMU directions   
        n_Cartesian = -vSW_Cartesian_KSMU;    
        Bfield_Cartesian = B_Cartesian_KSMU;  
        vSW_Cartesian = vSW_Cartesian_KSMU;

    else
        
    % - Conversion: KSM-U -> StPoint Frame
        SSN = ParaSystem.SSN;
        SSN_Direction = SSN.TiltedDipole.KSM_Rp.' ./ norm (SSN.TiltedDipole.KSM_Rp);
        ex_SSN = SSN_Direction;

        MAligned = ParaSystem.M .* [0; 0; 1]; 
        MTilted = ParaSystem.Tilt.RotationPhi*( ParaSystem.Tilt.RotationTheta*MAligned);
        ey_SSN = cross(MTilted, ex_SSN) ./ norm(cross(MTilted, ex_SSN));

        ez_SSN = cross(ex_SSN, ey_SSN);
        M_StPoint2KSMU = [ex_SSN, ey_SSN, ez_SSN];

        M_KSMU2StPoint =  (M_StPoint2KSMU)^-1;
            
        % Bfield and Vsw: in StPoint frame
        B_Cartesian_StPoint = M_KSMU2StPoint*B_Cartesian_KSMU;
        vSW_Cartesian_StPoint = M_KSMU2StPoint*vSW_Cartesian_KSMU;
      
        if Theta == 0
                    
            % In KSMU directions   
            n_Cartesian = -vSW_Cartesian_KSMU;    
            Bfield_Cartesian = B_Cartesian_KSMU;  
            vSW_Cartesian = vSW_Cartesian_KSMU;
    
        else

                    P_cart2sph = [ sin(Theta)*cos(Phi)     sin(Theta)*sin(Phi)    cos(Theta)    ;...
                                    cos(Theta)*cos(Phi)     cos(Theta)*sin(Phi)    -sin(Theta)   ;...
                                    -sin(Phi)               cos(Phi)               0             ...
                                ];
                    n_vec = [ 1 ; (-1/r)*drdtheta ; -1/(r*(sin(Theta)))*drdphi ].';
                    n_Cartesian_MB = (P_cart2sph.') * n_vec.' ;
                    n_Cartesian_StPoint = [n_Cartesian_MB(3), n_Cartesian_MB(1), n_Cartesian_MB(2)];
                    n_Cartesian_StPoint =  n_Cartesian_StPoint ./ norm(n_Cartesian_StPoint);
%                     n_Cartesian_StPoint = M_KSMU2StPoint*n_Cartesian_StPoint.';
                
                    % In SSN-centered directions   
                    n_Cartesian = n_Cartesian_StPoint;
                    Bfield_Cartesian = B_Cartesian_StPoint;
                    vSW_Cartesian = vSW_Cartesian_StPoint;
        end
    end
     
    PressureMag =  norm(cross(n_Cartesian, Bfield_Cartesian))*sqrt(1+Beta);
    PressureSW = -(1/2)*dot(n_Cartesian, vSW_Cartesian);
    PressureDifference = PressureMag - PressureSW;

    PressureAverage = (1/2)*(PressureMag+PressureSW);

    % PressureBalance = PressureDifference / (PressureAverage);
    PressureBalance = (PressureDifference);

end