function BField_Spectral = BField_Spectral(CoordinateDetails, PointCoordinates, ParaSystem)
    
    if strcmp(CoordinateDetails.Units, 'Rp')
        ConvFactor = ParaSystem.Rp ./ ParaSystem.r0;
    elseif strcmp(CoordinateDetails.Units, 'r0')
        ConvFactor = 1;
    end

%--- Extracting KSM Cartesian coordinates
    if strcmp(CoordinateDetails.Type, 'Cartesian_KSM')

        PositionVector = PointCoordinates .* ConvFactor;
        XKSM = PositionVector(1, :);
        YKSM = PositionVector(2, :);
        ZKSM = PositionVector(3, :);
        
    elseif strcmp(CoordinateDetails.Type, 'Spherical_MB')
           
        r = PointCoordinates(:, 3);
        Phi = PointCoordinates(:, 1);
        Theta = PointCoordinates(:, 2);

        ZMB = r .* cos(Theta);
        XMB = r .* sin(Theta) .* cos(Phi);
        YMB = r .* sin(Theta) .* sin(Phi);
        
        XKSM = ZMB;
        YKSM = XMB;
        ZKSM = YMB;
        
        PositionVector = [XKSM, YKSM, ZKSM] .* ConvFactor;
            
    end

%--- Position Vectors: in KSM-U Frame
    r = sqrt(XKSM.^2 + YKSM.^2 + ZKSM);
    er = PositionVector ./ r;

%--- Dipole Tilt: Rotation Matrices
    RotationPhi = ParaSystem.Tilt.RotationPhi;
    RotationTheta = ParaSystem.Tilt.RotationTheta;

%--- Expression for Tilted Dipole
    MAligned = ParaSystem.M .* [0; 0; -1]; 
    MTilted = RotationPhi*(RotationTheta*MAligned);
    
    MTiltedVector = ones(size(PositionVector));
    MTiltedVector(:, 1) = MTilted(1) .* MTiltedVector(:, 1);
    MTiltedVector(:, 2) = MTilted(2) .* MTiltedVector(:, 2);
    MTiltedVector(:, 3) = MTilted(3) .* MTiltedVector(:, 3);
    
    erKSM = er;
    M_Dot_er = MAligned(1).*erKSM(:, 1) + MAligned(2).*erKSM(:, 2) + MAligned(3).*erKSM(:, 3);
    BField_Spectral = (1./r.^3) .* ( 3*M_Dot_er .* erKSM -  MTiltedVector); 
     
%     Factor = ParaSystem.Rp ./ ParaSystem.r0;
%     PositionVector = PositionVector ./ Factor;
%     scatter3(PositionVector(1), PositionVector(2), PositionVector(3), 80, 'filled')
%     
%     PointCoordinates = PointCoordinates ./ Factor;
%     scatter3(PointCoordinates(1), PointCoordinates(2), PointCoordinates(3), 80, 'filled')
%     
          
end