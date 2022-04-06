
function ParaSystem2 = StagnationPoint(ParaSystem)

AlphaDipole = ParaSystem.AlphaDipole;
ParaSystem2 = ParaSystem;

% ----  Finding B.v=0 locus in NMM meridian plane

RotationPhi = ParaSystem.Tilt.RotationPhi;
RotationTheta = ParaSystem.Tilt.RotationTheta;
RotationAlphaDipole = ParaSystem.Tilt.RotationAlphaDipole;            

    CoordinateDetails.Units = 'Rp';
    CoordinateDetails.Type = 'Cartesian_KSM';
    CoordinateDetails.Component = 'BdotV';
    XGuess_List = (15:1:35);
    XYZ_SolutionList = ones(3, length(XGuess_List));
    for kX = 1:length(XGuess_List)
        
        XGuess = XGuess_List(kX);
        XYZ_Guess = [XGuess, 0, 0];
        if abs(ParaSystem.Tilt.Phi)<pi/2
            XYZ_Guess = (RotationTheta*XYZ_Guess.');
        else
            XYZ_Guess = (RotationTheta.'*XYZ_Guess.');
        end
        XYZ_Guess = (RotationAlphaDipole.'*XYZ_Guess);
        
        NormalField = @(Z) BdotV(CoordinateDetails, horzcat([XYZ_Guess(1); Z*tan(AlphaDipole); Z]), ParaSystem);
        
        opts = optimoptions('fsolve', 'algorithm','levenberg-marquardt');
        Z_Sol = fsolve(NormalField, XYZ_Guess(3), opts);
        XYZ_Sol = [XYZ_Guess(1), Z_Sol* tan(AlphaDipole), Z_Sol];
        BdotV(CoordinateDetails, horzcat([XYZ_Sol(1); XYZ_Sol(2); XYZ_Sol(3)]), ParaSystem)
        XYZ_SolutionList(:,kX) = XYZ_Sol;

    end

    Interpolant_NormalField_XToZ = griddedInterpolant(XYZ_SolutionList(1,:), XYZ_SolutionList(3,:), 'spline');
    Interpolant_NormalField_XToY = griddedInterpolant(XYZ_SolutionList(1,:), XYZ_SolutionList(2,:), 'spline');

    BDotV = plot3(XYZ_SolutionList(1,:), XYZ_SolutionList(2,:), XYZ_SolutionList(3,:), 'linewidth', 3);
    Colour_Equator = cbrewer2('Greens', 10);
    set(BDotV, 'Color', horzcat(Colour_Equator(8,:), 0.3), 'LineStyle', '-.')


% ----  Finding SSN on B.v=0 locus in NMM meridian plane

    drdtheta = 0;
    drdphi = 0;
    XGuess = 20;
    SSN_Guess = [XGuess, Interpolant_NormalField_XToY(XGuess),  Interpolant_NormalField_XToZ(XGuess)];

    CoordinateDetails.SSN = 1;
    funPB_SNN = @(X) PressureBalance_Equation(CoordinateDetails, horzcat([X; Interpolant_NormalField_XToY(X); Interpolant_NormalField_XToZ(X) ]), drdtheta, drdphi, ParaSystem);

    opts = optimoptions('fsolve', 'algorithm', 'levenberg-marquardt', 'FunctionTolerance', 1.e-10);
    SSN_Solution_X = fsolve( funPB_SNN, SSN_Guess(1), opts);

    SSN.TiltedDipole.KSM_Rp = [SSN_Solution_X, Interpolant_NormalField_XToY(SSN_Solution_X), Interpolant_NormalField_XToZ(SSN_Solution_X)];

    funPB_SNN(SSN.TiltedDipole.KSM_Rp(1))
    BdotV(CoordinateDetails, horzcat([SSN.TiltedDipole.KSM_Rp(1); SSN.TiltedDipole.KSM_Rp(2); SSN.TiltedDipole.KSM_Rp(3)]), ParaSystem)

    Radius = 0.4;
    [u,v,w]=sphere(25);
    SSNPlot = surf(SSN.TiltedDipole.KSM_Rp(1)+Radius*u,SSN.TiltedDipole.KSM_Rp(2)+Radius*v,SSN.TiltedDipole.KSM_Rp(3)+Radius*w);
    set(SSNPlot, 'FaceLighting', 'none', 'FaceColor', 'red', 'EdgeColor', 'none', 'FaceAlpha', 0.5)

% ----  Translating magnetic equatorial plane into SSN-centered phi=0 plane
%     Equ_Plane.Vertices = Equ_Plane.Vertices + [SSN.TiltedDipole.KSM_Rp; SSN.TiltedDipole.KSM_Rp; SSN.TiltedDipole.KSM_Rp; SSN.TiltedDipole.KSM_Rp];

    ParaSystem2.SSN = SSN;

end
