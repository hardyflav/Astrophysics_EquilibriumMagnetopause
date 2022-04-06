function BField = BField(CoordinateDetails, PointCoordinates, ParaSystem)
    
    if strcmp(CoordinateDetails.Units, 'Rp')
        ConvFactor = ParaSystem.Rp ./ ParaSystem.r0;
    elseif strcmp(CoordinateDetails.Units, 'r0')
        ConvFactor = 1;
    end

%--- Extracting KSM Cartesian coordinates
    if strcmp(CoordinateDetails.Type, 'Cartesian_KSM')

        PositionVector = PointCoordinates .* ConvFactor;
        
    elseif strcmp(CoordinateDetails.Type, 'Spherical_MB')
           
        r = PointCoordinates(3, 1);
        Phi = PointCoordinates(1, 1);
        Theta = PointCoordinates(2, 1);

        ZMB = r * cos(Theta);
        XMB = r * sin(Theta) * cos(Phi);
        YMB = r * sin(Theta) * sin(Phi);
        
        XKSM = ZMB;
        YKSM = XMB;
        ZKSM = YMB;
        
        PositionVector = [XKSM; YKSM; ZKSM] .* ConvFactor;
            
    end

%--- Position Vectors: in KSM-U Frame
    r = norm(PositionVector);
    er = PositionVector ./ r;

%--- Dipole Tilt: Rotation Matrices
    RotationPhi = ParaSystem.Tilt.RotationPhi;
    RotationTheta = ParaSystem.Tilt.RotationTheta;

%--- Expression for Tilted Dipole
    MAligned = ParaSystem.M .* [0; 0; -1]; 
    MTilted = RotationPhi*(RotationTheta*MAligned);
    erKSM = er;
    BField = (1/r^3) .* ( 3*dot(MTilted, erKSM) .* erKSM -  MTilted); 
     
%     Factor = ParaSystem.Rp ./ ParaSystem.r0;
%     PositionVector = PositionVector ./ Factor;
%     scatter3(PositionVector(1), PositionVector(2), PositionVector(3), 80, 'filled')
%     
%     PointCoordinates = PointCoordinates ./ Factor;
%     scatter3(PointCoordinates(1), PointCoordinates(2), PointCoordinates(3), 80, 'filled')
%     
          
end