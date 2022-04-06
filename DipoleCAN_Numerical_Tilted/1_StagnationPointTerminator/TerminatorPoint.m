function ParaSystem2 = TerminatorPoint(ParaSystem)

    AlphaDipole = ParaSystem.AlphaDipole;
    ParaSystem2 = ParaSystem;

% ----  Finding A2 on norm(B x v) = 0 locus in NMM meridian plane
    CoordinateDetails.Units = 'Rp';
    CoordinateDetails.Type = 'Cartesian_KSM';
    CoordinateDetails.Component = 'BdotV';

%     ZGuess_List = (15:1:35);
    ZGuess_List = (-35:1:-15);
    XYZ_SolutionList = ones(3, length(ZGuess_List));
    for kZ = 1:length(ZGuess_List)
        
        ZGuess = ZGuess_List(kZ);
        XYZ_Guess = [0, ZGuess*tan(AlphaDipole), ZGuess];
    
        TangentialField = @(X) BdotV(CoordinateDetails, horzcat([X; XYZ_Guess(3)*tan(AlphaDipole); XYZ_Guess(3)]), ParaSystem);
        opts = optimoptions('fsolve', 'algorithm','levenberg-marquardt');
        X_Sol = fsolve(TangentialField, XYZ_Guess(1), opts);
        XYZ_Sol = [X_Sol, XYZ_Guess(2), XYZ_Guess(3)];
        
        BdotV(CoordinateDetails, horzcat([XYZ_Sol(1); XYZ_Sol(2); XYZ_Sol(3)]), ParaSystem)
        XYZ_SolutionList(:,kZ) = XYZ_Sol;
        
    end
    
    Interpolant_NormalField2_ZToX = griddedInterpolant(XYZ_SolutionList(3,:), XYZ_SolutionList(1,:), 'spline');
    Interpolant_NormalField2_ZtoY = griddedInterpolant(XYZ_SolutionList(3,:), XYZ_SolutionList(2,:), 'spline');

    BDotV2 = plot3(XYZ_SolutionList(1,:), XYZ_SolutionList(2,:), XYZ_SolutionList(3,:), 'linewidth', 3);
    Colour_Equator = cbrewer2('Greens', 10);
    set(BDotV2, 'Color', horzcat(Colour_Equator(8,:), 0.3), 'LineStyle', '-.')
        
    
% ----  Finding A2 on (Bxv)=0 locus in NMM meridian plane

%     ZGuess = 10;
    ZGuess = -10;
    A1_Guess = [Interpolant_NormalField2_ZToX(ZGuess), Interpolant_NormalField2_ZtoY(ZGuess),  ZGuess];

    CoordinateDetails.Units = 'Rp';
    CoordinateDetails.Type = 'Cartesian_KSM';
    funA2 = @(Z) norm(BField(CoordinateDetails, [Interpolant_NormalField2_ZToX(Z); Interpolant_NormalField2_ZtoY(Z); Z], ParaSystem)) - 1/2;    

    funA2(ZGuess);
    BdotV(CoordinateDetails, horzcat([A1_Guess(1); A1_Guess(2); A1_Guess(3)]), ParaSystem)

    opts = optimoptions('fsolve', 'algorithm', 'levenberg-marquardt', 'FunctionTolerance', 1.e-10);
    A1_Solution_Z = fsolve( funA2, A1_Guess(3), opts);
    A1_Solution = [Interpolant_NormalField2_ZToX(A1_Solution_Z), Interpolant_NormalField2_ZtoY(A1_Solution_Z), A1_Solution_Z];
    funA2(A1_Solution(3));
    BdotV(CoordinateDetails, horzcat([A1_Solution(1); A1_Solution(2); A1_Solution(3)]), ParaSystem);

    Radius = 0.4;
    [u,v,w]=sphere(25);
    A1Plot = surf(A1_Solution(1)+Radius*u,A1_Solution(2)+Radius*v,A1_Solution(3)+Radius*w);
    set(A1Plot, 'FaceLighting', 'none', 'FaceColor', 'red', 'EdgeColor', 'none', 'FaceAlpha', 0.5)

    norm(A1_Solution);

    
A1.TiltedDipole.KSM_Rp = A1_Solution;
ParaSystem2.A1 = A1;

end