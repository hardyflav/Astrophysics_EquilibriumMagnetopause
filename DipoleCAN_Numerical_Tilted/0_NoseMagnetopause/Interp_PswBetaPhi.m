clear variables
clc;
p = genpath('/Users/flavienhardy/Documents/Git/MagnetopauseTilt');
addpath(p);

global DegToRad
global RadToDeg
DegToRad = pi/180;
RadToDeg = 180/pi;


%% Setting up system
Inclination = 1;
Planet = 'Saturn';
ParaSystem = System_ini(Planet, Inclination);
ParaSystem.Inclination = Inclination;



%% Computing CAN-Disk field


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

% ----  Packaging Output
    BCAN = BCAN_Rotated;
    ParaSystem.BCAN = BCAN;


%% Setting up grid
PswList = (0.001: 0.001 : 10^(-0.5)).';
RmpList = (15:1:35);
PHIList = (0:1:180).' .* DegToRad;

[PHIGrid, RmpGrid, PswGrid] = ndgrid(PHIList, RmpList, PswList);
BetaList = 0.* PswGrid;

for kPhi = 1:length(PHIList)
    ParaSystem.Phi = PHIList(kPhi);
    
    for krMP = 1:length(RmpList)
        ParaSystem.krMP = RmpList(krMP);
        ParaSystem.a = 0.030 * ParaSystem.krMP + 6.1; % Rs
        ParaSystem.b = 0.63 * ParaSystem.krMP + 4.8; % Rs
        ParaSystem.mu0I = 1.12 * ParaSystem.krMP + 26.7; % nT
        
        
        for kPsw = 1:length(PswList)
            ParaSystem.Psw = PswList(kPsw);
            
            ParaSystem_InterpPswBetaPhi = System_InterpPswBetaPhi(ParaSystem);            
            FuncBeta = @(X) PB_Nose_BetaDetermination(ParaSystem.krMP, X(1), X(2), ParaSystem_InterpPswBetaPhi);
            
            AlphaList = (-30:0.1:30).' * pi/180;
            for kAlpha = 1:length(AlphaList)

                Alpha = AlphaList(kAlpha);
                PBFunctionR = @(Beta) FuncBeta([Alpha, Beta]);
                rSol = fsolve(PBFunctionR, 1);
                BetaList(kPhi, krMP, kPsw) = rSol;

            end
            
            
            
        end
    end
    
end



