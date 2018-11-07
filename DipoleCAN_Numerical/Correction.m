function [SurfaceCorrected] = Correction(rSubSolarNose, DeltaThetaDeg, DeltaPhiDeg, ThetaCusp, ThetaMaxDeg, PhiMaxDeg, rEquator, rMeridian, SurfaceTot, NumIterationsTop, NumIterationsBottom, SystemParameters)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs
%    - DeltaThetaDeg, DeltaPhiDeg: scalars, angular increments, in degrees
%    - ThetaCusp: scalar, theta value of the cusp, in degrees
%    - ThetaMaxDeg, PhiMaxDeg: scalars, maximum values for theta and phi,
%           in degrees
%    - rEquator: vector, solution in the equatorial plane
%           theta = 0:DeltaThetaDeg:ThetaMaxDeg
%    - rMeridian: vector, solution in the noon-midnight meridian plane
%           theta = 0:DeltaThetaDeg:ThetaMaxDeg
%    - SurfaceTot: phi*theta array, surface to be corrected
%           on the entire grid
%    - NumIterationsTop, NumIterationsBottom: scalars, number of iterations
%           for the correction of the top/bottom sub-grid
%
% Outputs
%     - SurfaceCorrected: phi*theta array, corrected surface
%           on the entire grid
%     - rCroppedCorrected: vector, corrected upper sub-grid
%     - SurfaceCroppedExtrapolated: phi*theta array, corrected surface on
%     the sub-upper grid
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% Cropping BEFORE cusp
    ThetaMaxDegSubSolar =  (ThetaCusp -  mod(ThetaCusp, DeltaThetaDeg)) - 0*DeltaThetaDeg;
    PhiMaxDegSubSolar = PhiMaxDeg;

    [~, NbPointsGridSubSolarInside, NbPointsThetaSubSolarInside, NbPointsPhiSubSolarInside, ~, ~] = GridDetails(ThetaMaxDegSubSolar, PhiMaxDegSubSolar, DeltaThetaDeg, DeltaPhiDeg);

    ArrayGridInside = SurfaceTot(2:end-1, 2:end-1);
    ArrayGridTot = SurfaceTot;
    ArrayGridSubSolarInside = ArrayGridInside(1:end, 1:NbPointsThetaSubSolarInside);

    % theta = DeltaThetaDeg : DeltaThetaDeg : ThetaMaxDegSubSolar
    % Phi = DeltaPhiDeg : DeltaPhiDeg : PhiMaxDegSubSolar
    % Size = NbPointsGridSubSolarInside = (NbPointsThetaSubSolar-2) * (NbPointsPhiSubSolar-2)
    rGuess = ArrayGridSubSolarInside(:);


%% Correcting sub-solar cap
    epsilon = 1.e-10;
    NumIterations = 15;

    opts = optimoptions('fsolve','algorithm','Levenberg-Marquardt', 'Display', 'iter');
    opts.MaxIterations = NumIterations;
    opts.MaxFunEvals = 100;
    opts.SpecifyObjectiveGradient = true;
    opts.CheckGradients = false;
    opts.FiniteDifferenceType = 'central';
    [rNew] = fsolve( @(r) F_valuesJacobian_Numerical(r, rEquator, rMeridian, rSubSolarNose, ThetaMaxDegSubSolar, PhiMaxDegSubSolar, DeltaThetaDeg, DeltaPhiDeg, NbPointsGridSubSolarInside, SystemParameters, epsilon), rGuess, opts);
    rSubSolarCorrected = rNew;


%% Inclusion of the Meridian and Equator
    SubSolarGridCorrected = reshape(rNew,[NbPointsPhiSubSolarInside, NbPointsThetaSubSolarInside]);
    EquatorCrop = rEquator(2:NbPointsThetaSubSolarInside+1).';
    MeridianCrop = rMeridian(2:NbPointsThetaSubSolarInside+1).';

    SubSolarGridCorrectedCompletedTemp = vertcat(EquatorCrop, SubSolarGridCorrected, MeridianCrop);
    SubSolarNose = rSubSolarNose*ones(NbPointsPhiSubSolarInside+2, 1);
    GridCropCompleteCorrected = horzcat(SubSolarNose, SubSolarGridCorrectedCompletedTemp);



%% Manual linear extrapolation until cusp
    Offset = (ThetaCusp - ThetaMaxDegSubSolar + DeltaThetaDeg)/DeltaThetaDeg;
    ExtraRows = NaN*ones(NbPointsPhiSubSolarInside+2, Offset);
    SubSolarCapLong = horzcat(GridCropCompleteCorrected, ExtraRows);
    SubSolarCapLong(end, 1:size(SubSolarCapLong, 2)) = rMeridian(1:size(SubSolarCapLong, 2));
    SubSolarCapLong(1, 1:size(SubSolarCapLong, 2)) = rEquator(1:size(SubSolarCapLong, 2));

    for k = 1:Offset
        SubSolarCapLong(2:end-1, size(GridCropCompleteCorrected,2)+k) = 2*SubSolarCapLong(2:end-1, size(GridCropCompleteCorrected,2)+k-1) - SubSolarCapLong(2:end-1, size(GridCropCompleteCorrected,2)+k-2);
    end

    SurfaceCroppedExtrapolated = SubSolarCapLong;

    
%% Correcting the bottom grid: building the initial guess
% New left boundary
    % theta = 0:ThetaMax
    % phi = 0:PhiMax
    LeftBoundaryVector = SubSolarCapLong(:, end);

% New interior grid
    % theta = SubSolarCapLong(end)+delta : ThetaMax-delta (pour delta = 1
    % deg, theta = 71:124
    % phi = delta:max-delta
    BottomGrid = ArrayGridTot(2:end-1, size(SubSolarCapLong, 2)+1:end-1);

    TransitionWidth = (15 - mod(15,DeltaThetaDeg))/DeltaThetaDeg;
    TransitionWidthDeg = TransitionWidth*DeltaThetaDeg;
    GridRight = ArrayGridTot(:, size(SubSolarCapLong, 2)+TransitionWidth:end);

    ThetaInterpDeg_Dummy = ThetaCusp+TransitionWidthDeg : DeltaThetaDeg : ThetaMaxDeg;
    ThetaInterpDeg_Dummy2 = [ThetaCusp, ThetaInterpDeg_Dummy];
    ThetaInterp = ThetaInterpDeg_Dummy2*pi/180;

    PhiInterp = (0:DeltaPhiDeg:PhiMaxDeg)*pi/180;

    [ThetaInterp, PhiInterp] = ndgrid(ThetaInterp, PhiInterp);

    ValuesInterp = horzcat(LeftBoundaryVector, GridRight);

    F = griddedInterpolant(ThetaInterp, PhiInterp, ValuesInterp.', 'spline');

    Thetaq_Temp = (ThetaCusp:DeltaThetaDeg:ThetaMaxDeg)*pi/180;
    Phiq_Temp = (0:DeltaPhiDeg:PhiMaxDeg)*pi/180;
    [Thetaq, Phiq] = ndgrid(Thetaq_Temp, Phiq_Temp);
    Vq = F(Thetaq, Phiq).';

    InterpolatedGridInside = Vq((2:end-1), (2:end-1));
    InterpolatedGridTot = Vq((1:end), (1:end));
    rGuessBottom = InterpolatedGridInside(:);


%% Correcting the bottom grid
    NbPointsGridBottomInside = size(rGuessBottom,1) ;
    NbPointsPhiBottomInside = size(BottomGrid, 1);
    NbPointsThetaBottomInside = size(BottomGrid, 2);

    NumIterationsBottom = 15;
    opts = optimoptions('fsolve','algorithm','Levenberg-Marquardt', 'Display', 'iter');
    opts.MaxIterations = NumIterationsBottom;
    opts.MaxFunEvals = 100000;
    opts.SpecifyObjectiveGradient = true;
    opts.CheckGradients = false;
    opts.FiniteDifferenceType = 'central';

    rLeft = LeftBoundaryVector;
    rBottom = rEquator( (ThetaCusp+DeltaThetaDeg)/DeltaThetaDeg   : end  );
    rTop = rMeridian( (ThetaCusp+DeltaThetaDeg)/DeltaThetaDeg   : end  ).';

    rNew = fsolve( @(r) Values_Jacobian_Bottom( r, rBottom, rLeft, rTop, ThetaCusp, ThetaMaxDeg, PhiMaxDeg, DeltaThetaDeg, DeltaPhiDeg, NbPointsGridBottomInside, SystemParameters, epsilon), rGuessBottom, opts);
    rTailRegionCorrected = rNew;


%% Merging the two sub-surfaces
    BottomGridCorrected_temp1 = reshape(rTailRegionCorrected, [NbPointsPhiBottomInside, NbPointsThetaBottomInside]);
    BottomGridCorrected_temp2 = vertcat(rBottom(2:end-1).', BottomGridCorrected_temp1, rTop(2:end-1));
    BottomGridCorrected_temp3 = horzcat(rLeft, BottomGridCorrected_temp2);

    rRightCorrected = 2*BottomGridCorrected_temp3(:,end) - BottomGridCorrected_temp3(:,end-1);

    BottomGridCorrected = horzcat(BottomGridCorrected_temp3, rRightCorrected);
    % theta = rLeftBoundary : delta : ThetaMax (for delta = 1deg, 70:125)
    % phi = 0:90
    SurfaceCorrected = horzcat(SubSolarCapLong, BottomGridCorrected(:,2:end));



end
