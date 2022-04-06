%{
function Psw = PhiRmp2Psw( PhiTilt, rMP )
    

%% Building Array: PhiTilt - Psw - rMP

%--- rMP range at Saturn
    Pswnpa = 0.2;
    PhiList = [0:20:180]; 
    rMPList = (15:2:35);
    PswEff_Array = [...
        [ 0.2/(1+0.584);  0.2/(1+1.08); 0.1/(1+1.8276)], ...
                ] ;
%}            
            
function rMPDiff = testPsw(PhiTilt, rMP, Psw)

    Inclination = 1;
    Planet = 'Saturn';
    Conditions_PhiTilt_Psw_rMP = [ PhiTilt, Psw, rMP  ];
    ParaSystem = System_ini(Planet, Inclination, Conditions_PhiTilt_Psw_rMP);
    ParaSystem.BShielding.X = 0;
    ParaSystem.BShielding.Y = 0;
    ParaSystem.BShielding.Z = 0;

    
%% ------------------------------------------------------------------------
% Computing ROTATED CAN Disk X-Y-Z KSM Components FROM KSM coordinates

% ----  Setting up Cartesian grid for CAN Disk Field Interpolation
    ParaSystem.OnOff_CAN = 1;
    CAN_DiskParameters = ParaSystem.CAN_DiskParameters;

    XList_KSM = (-36 : 2 : 36) .';   
    YList_KSM = (-36 : 2 : 36) .';   
    ZList_KSM = (-36 : 2 : 36) .';   
    [Position_KSM] = [XList_KSM, YList_KSM, ZList_KSM] .';

    [XGrid_KSM, YGrid_KSM, ZGrid_KSM] = ndgrid(XList_KSM, YList_KSM, ZList_KSM);

% ----  Rotating back to dipole-aligned case: Position_Aligned = [X; Y; Z]
    P_AlignedToRot = ParaSystem.Tilt.RotationPhi * ParaSystem.Tilt.RotationTheta;
    P_RotToAligned = P_AlignedToRot.';

    XGrid_Aligned = P_RotToAligned(1,1) .* XGrid_KSM + P_RotToAligned(1,2) .* YGrid_KSM  + P_RotToAligned(1,3) .* ZGrid_KSM;
    YGrid_Aligned = P_RotToAligned(2,1) .* XGrid_KSM + P_RotToAligned(2,2) .* YGrid_KSM  + P_RotToAligned(2,3) .* ZGrid_KSM;
    ZGrid_Aligned = P_RotToAligned(3,1) .* XGrid_KSM + P_RotToAligned(3,2) .* YGrid_KSM  + P_RotToAligned(3,3) .* ZGrid_KSM;

    ZGrid_Aligned_MB = XGrid_Aligned;
    XGrid_Aligned_MB = YGrid_Aligned;
    YGrid_Aligned_MB = ZGrid_Aligned;

    [Phi_Aligned, Elevation_Aligned, r_Aligned] = cart2sph( XGrid_Aligned_MB, YGrid_Aligned_MB, ZGrid_Aligned_MB );
    Theta_Aligned = pi/2 - Elevation_Aligned;
    Phi_Aligned = wrapToPi(Phi_Aligned);

    BCAN_Aligned = B_CAN_Corrected(CAN_DiskParameters, r_Aligned, Theta_Aligned, Phi_Aligned);
    
% ----  Rotating to initial tilted-dipole case: Position_Tilted = [XKSM; YKSM; ZKSM]
    BCAN_Rotated.X = P_AlignedToRot(1,1) .* BCAN_Aligned.X + P_AlignedToRot(1,2) .* BCAN_Aligned.Y  + P_AlignedToRot(1,3) .* BCAN_Aligned.Z;
    BCAN_Rotated.Y = P_AlignedToRot(2,1) .* BCAN_Aligned.X + P_AlignedToRot(2,2) .* BCAN_Aligned.Y  + P_AlignedToRot(2,3) .* BCAN_Aligned.Z;
    BCAN_Rotated.Z = P_AlignedToRot(3,1) .* BCAN_Aligned.X + P_AlignedToRot(3,2) .* BCAN_Aligned.Y  + P_AlignedToRot(3,3) .* BCAN_Aligned.Z;

    BCAN_Rotated.XInterp = griddedInterpolant(XGrid_KSM, YGrid_KSM, ZGrid_KSM, BCAN_Rotated.X, 'spline');
    BCAN_Rotated.YInterp = griddedInterpolant(XGrid_KSM, YGrid_KSM, ZGrid_KSM, BCAN_Rotated.Y, 'spline');
    BCAN_Rotated.ZInterp = griddedInterpolant(XGrid_KSM, YGrid_KSM, ZGrid_KSM, BCAN_Rotated.Z, 'spline');

% ----  Optional: Visualisation of CAN-disk field 
%{
    XList = (-20:1:30);
    YList = (-20:1:20);
    ZList = (-20:1:20);

    [X, Y, Z] = ndgrid( XList, YList, ZList );
    figure;
        PlotQuiver = quiver3(X, Y, Z, BCAN_Rotated.XInterp(X, Y, Z), BCAN_Rotated.YInterp(X, Y, Z), BCAN_Rotated.ZInterp(X, Y, Z));
        set(PlotQuiver, 'AutoScaleFactor', 2, 'LineWidth', 1);
        axis equal
%}

% ----  Packaging Output
    BCAN = BCAN_Rotated;
    ParaSystem.BCAN = BCAN;
    
    
    
%% -------------------------------------------------------------------------------
% 3D Plot of the System

    RotationOptions = ["DipoleTilted", "DipoleAligned"];
    Rotation = RotationOptions(1);
%     Plot = Plot3D_ini(Rotation, ParaSystem);

    
    
%% -------------------------------------------------------------------------------
% Finding Position of Anchor Points: Stagnation Point and Terminator Point

% ----  Finding B.v=0 locus in NMM meridian plane + Stagnation Point
    ParaSystem = NoseDetermination(ParaSystem);

% ----  Finding Bxv=0 locus in NMM meridian plane + Terminator Point
    ParaSystem = TerminatorPoint(ParaSystem);

    rMPDiff = rMP - norm(ParaSystem.Nose.KSMU_Rp);
                
end            
            
            
            
            
            
            
 %{           
            
            
            

    krMP = 1;        
    NoseFunction = @(X) PB_Nose_BetaDetermination(rMPList(krMP), X(1), X(2), ParaSystem);

    AlphaList = (-30:0.1:30).' * pi/180;
    BetaList = 0 * AlphaList;
    
%         kAlpha=40
    for kAlpha = 1:length(AlphaList)
        
        Alpha = AlphaList(kAlpha);
        PBFunctionR = @(Beta) NoseFunction([Alpha, Beta]);
        rSol = fsolve(PBFunctionR, 1);
        BetaList(kAlpha) = rSol;
        
    end            
            
            
            
            
            
            
            
            
            
            
            
            
            
  %{          

%--- Solar Wind pressure range at Saturn
    PswList = 10.^( -3:0.5:-1 ) .';
    PhiList = [0; 20; 40; 60; 80; 90; 110; 120; 130; 140; 160; 180]; 

    rMP_Array = [...
        [ 33.32; 27.61; 23.04; 19.2945; 16.491 ],...
        [ 33.09; 27.4219; 22.8924; 19.1744; 16.4621 ], ...
        [ 32.7753; 27.1787; 22.7120; 19.0127; 16.5187 ], ...
        [ 32.5495; 27.0127; 22.6034; 18.9582; 15.9145 ], ...
        [ 32.5078; 26.9959; 22.6138; 18.9945; 16.7287 ], ...
        [ 32.5603; 27.0472; 22.6596; 19.0451; 16.8207 ], ...
        [ 32.7747; 27.2215; 22.7922; 19.1362; 16.8839 ], ...
        [ 32.9137; 27.3309; 22.8751; 19.2256; 16.8483 ], ...
        [ 33.0557; 27.4393; 22.9479; 19.2604; 16.8256 ], ...
        [ 33.1866; 27.5377; 23.0178; 19.3247; 16.7286 ], ...
        [ 33.3562; 27.6507; 23.0847; 19.3550; 16.6404 ], ...
        [ 33.3227; 27.6132; 23.0389; 19.2945; 16.4910 ], ...
                ] ;

%--- Filling-in array with rMP values
    [PhiGrid, PswGrid] = ndgrid(PhiList, PswList);
    rMPGRid = 0.*PhiGrid;

    for kPhi = 1:length(PhiList)
        PhiVal = PhiList(kPhi);
        rMPVal = rMP_Array(:, kPhi);
        rMPGRid( PhiGrid==PhiVal ) = rMPVal;
    end

%--- Using array to determine SPECIFIC rMP-contour and interpolant
    Contours_rMP = contourc( PhiList, log10(PswList), rMPGRid.', [rMP, rMP] );
    
%--- Finding SPECIFIC Psw value for chosen (PhiTilt, rMP)
    rMP_Interp = griddedInterpolant( Contours_rMP(1, 2:end),  Contours_rMP(2, 2:end));
    Psw = 10^rMP_Interp(PhiTilt);
    
   %}
                
                
end
%}


    