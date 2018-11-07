clear all;
clc;
p = genpath('/Users/flavien/Documents/Github/EquilibriumMagnetopause/DipoleCAN_Numerical');
addpath(p)


%% Environment Definition: Planet and Pressure Balance Parameters
    Planet = "Saturn";
    SystemParameters = SystemParameters(Planet);


%% Grid Definition for Surface Optimisations
    ThetaMaxDeg = 120;                             % Maximum value of theta on the grid, degrees
    PhiMaxDeg = 90;                                % Maximum value of phi on the grid, degrees
    DeltaThetaDeg = 2;
    DeltaPhiDeg = 2;
    ThetaSpan = (0 : DeltaThetaDeg : (ThetaMaxDeg))*pi/180;
    PhiSpan = (0 : DeltaPhiDeg: PhiMaxDeg).*pi/180;
    [NbPointsGrid, NbPointsGridInside,      ...
    NbPointsThetaInside, NbPointsPhiInside, ...
    ThetaSpanVector, PhiSpanVector]         = GridDetails(ThetaMaxDeg, PhiMaxDeg, DeltaThetaDeg, DeltaPhiDeg);


%% Construction of the Initial Surface
    PlotsPossibilitiesInitial = ["CurveEquator", "CurveMeridianParts", "CurveMeridian", "SurfaceInitialTot", "SurfaceInitialTotFieldLines", "GridAnalysis", "GridWrapped", "GridWrappedFieldLines"];
    Plots =["CurveEquator", "CurveMeridian", "GridWrapped"];
    [rSubSolarNose, rEquatorInterpolant, rMeridianInterpolant, ThetaCuspCurves, InitialSurfaceInterpolant] = InitialGuess(SystemParameters, ThetaMaxDeg, PhiMaxDeg, DeltaThetaDeg, DeltaPhiDeg, Plots);

    rEquator = rEquatorInterpolant(ThetaSpan).';
    rMeridian = rMeridianInterpolant(ThetaSpan).';

    [PhiSpanOptimisation, ThetaSpanOptimisation] = ndgrid(PhiSpan, ThetaSpan);
    SurfaceInitialTot = InitialSurfaceInterpolant(PhiSpanOptimisation, ThetaSpanOptimisation);
    ThetaCusp = ThetaCuspCurves - mod(ThetaCuspCurves, DeltaThetaDeg);


%% Cropping and Correcting the Initial Surface
    NumIterationsTop = 15;
    NumIterationsBottom = 15;
    [SurfaceCorrected] = Correction(rSubSolarNose, DeltaThetaDeg, DeltaPhiDeg, ThetaCusp, ThetaMaxDeg, PhiMaxDeg, rEquator, rMeridian, SurfaceInitialTot, NumIterationsTop, NumIterationsBottom, SystemParameters);

    
%% Ploting the Results
    PlotsPossibilitiesCorrection= ["SurfaceCorrected", "SurfaceCorrectedFieldLines", "GridWrapped", "GridWrappedFieldLines"];
    PlotsCorrection = ["GridWrapped"];
    CorrectionPlot(PlotsCorrection, DeltaThetaDeg, DeltaPhiDeg, ThetaMaxDeg, PhiMaxDeg, rSubSolarNose, SurfaceCorrected, SystemParameters)

    
    
    
    
    
   

   
    
    
    
    
% DE CÔTÉ
%{
%% Inclusion of the Meridian and Equator
    SubSolarGridCorrected = reshape(rNew,[NbPointsPhiSubSolarInside, NbPointsThetaSubSolarInside]);
    EquatorCrop = rEquator(2:NbPointsThetaSubSolarInside+1).';
    MeridianCrop = rMeridian(2:NbPointsThetaSubSolarInside+1);

    SubSolarGridCorrectedCompletedTemp = vertcat(EquatorCrop, SubSolarGridCorrected, MeridianCrop.');
    SubSolarNose = rSubSolarNose*ones(NbPointsPhiSubSolarInside+2, 1);
    GridCropCompleteCorrected = horzcat(SubSolarNose, SubSolarGridCorrectedCompletedTemp);


%% Manual linear extrapolation until cusp
    Offset = (ThetaCusp - ThetaMaxDegSubSolar + DeltaThetaDeg)/DeltaThetaDeg;
    ExtraRows = NaN*ones(NbPointsPhiSubSolarInside+2, Offset);
    SubSolarCapLong = horzcat(GridCropCompleteCorrected, ExtraRows);
    SubSolarCapLong(end, 1:size(SubSolarCapLong, 2)) = rMeridian(1:size(SubSolarCapLong, 2)).';
    SubSolarCapLong(1, 1:size(SubSolarCapLong, 2)) = rEquator(1:size(SubSolarCapLong, 2)).';

    for k = 1:Offset
        SubSolarCapLong(2:end-1, size(GridCropCompleteCorrected,2)+k) = 2*SubSolarCapLong(2:end-1, size(GridCropCompleteCorrected,2)+k-1) - SubSolarCapLong(2:end-1, size(GridCropCompleteCorrected,2)+k-2);
    end

    SurfaceCroppedExtrapolated = SubSolarCapLong;



%% Plot corrected sub-solar cap

figure;

    hold on
    % Ring
    Bp = 20000*10^(-9);     % Equatorial field, T
    Rp = 60280*10^3;        % Planet radius, m
    mu_0 = 4*pi*10^(-7);
    Pswnpa = 0.02;
    b0 = sqrt(2*mu_0*Pswnpa*10^(-9));   % Field scale
    r0 = (2*Bp*Rp^3/b0)^(1/3);          % Distance scale
    RpScaled = Rp/r0;
    r = 66900*10^3/r0; %inner radius
    R = 140180*10^3/r0; %outer radius
    steps = 100;
    theta = 0:2*pi/steps:2*pi;
    r_outer = ones(size(theta))*R;
    r_inner = ones(size(theta))*r;
    rr = [r_outer;r_inner];
    theta = [theta;theta];
    xx = rr.*cos(theta);
    yy = rr.*sin(theta);
    zz = zeros(size(xx));
    Rings = surf(xx, yy,zz);
    set(Rings,'FaceColor',[0.6000    0.5000    0.5000],...
              'FaceAlpha', 0.5,...
              'LineStyle','none');

    % CAN Disk
    a = 8*Rp/r0;      % Inner radius of the current distribution (in r0)
    b = 15.5*Rp/r0;   % Outer radius of the current distribution (in r0)
    D = 3*Rp/r0;       % Disk Semi-thickness (in r0)
    t = linspace(0,2*pi);
    rin = a;
    rout = b;
    center = [0, 0];
    xin = rin*cos(t);
    xout = rout*cos(t);
    yin = rin*sin(t);
    yout = rout*sin(t);
    z1 = -D;
    z2 = D;
    
    bottom = patch(center(1)+[xout,xin], ...
               center(2)+[yout,yin], ...
               z1*ones(1,2*length(xout)),'', 'FaceAlpha', 0.5, 'LineStyle', 'none');
    top = patch(center(1)+[xout,xin], ...
            center(2)+[yout,yin], ...
            z2*ones(1,2*length(xout)),'', 'FaceAlpha', 0.5, 'LineStyle', 'none');
    [X,Y,Z] = cylinder(1,length(xin));
    outer = surf(rout*X+center(1), ...
             rout*Y+center(2), ...
             Z*(z2-z1)+z1, 'FaceAlpha', 0.4, 'LineStyle', 'none');
    inner = surf(rin*X+center(1), ...
             rin*Y+center(2), ...
             Z*(z2-z1)+z1, 'FaceAlpha', 0.4, 'LineStyle', 'none');

    
    ThetaSpanSubSolar = (0:DeltaThetaDeg:ThetaCusp)*pi/180;
    PhiSpanSubSolar = (0:DeltaPhiDeg:PhiMaxDegSubSolar)*pi/180;
    ConstructionSurfaceContours(SurfaceCroppedExtrapolated.', ThetaSpanSubSolar, PhiSpanSubSolar);
    set(gcf,'color','w');
    axis([-1.5, 1.5, -2.3, 2.3, -2.3, 2.3])
    box on 
    view([45+90+5 30])
    grid on
    
    
    
%% Analysis of the sub-solar cap
figure;
% SurfaceCroppedExtrapolatedBis = SurfaceCroppedExtrapolated(:);
% rSubSolarCap = SurfaceCroppedExtrapolatedBis(:);
% PB_Assessment = F_values_Numerical(rGuess, rEquator, rMeridian, rSubSolarNose, ThetaMaxDegSubSolar, PhiMaxDegSubSolar, DeltaThetaDeg, DeltaPhiDeg, NbPointsGridSubSolarInside, SystemParameters);
PB_Assessment = F_values_Numerical(rSubSolarCorrected, rEquator, rMeridian, rSubSolarNose, ThetaMaxDegSubSolar, PhiMaxDegSubSolar, DeltaThetaDeg, DeltaPhiDeg, NbPointsGridSubSolarInside, SystemParameters);
ErrorVector = log10(abs(PB_Assessment));
ErrorGrid = reshape(ErrorVector, length(PhiSpanSubSolar)-2, length(ThetaSpanSubSolar)-2);
pcolor(ErrorGrid);
colorbar




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

%%% TEST : 2D interpolation middle section
% CurvesDifference = rEquator-rMeridian;
% CurvesIntersectionDeg = 5 - mod(5, DeltaThetaDeg);
% 
% k = (CurvesIntersectionDeg)/DeltaThetaDeg-1;
% while CurvesDifference(k+2)*CurvesDifference(k+1)>0
%         CurvesIntersectionDeg = CurvesIntersectionDeg+DeltaThetaDeg;
%         k=k+1;
% end
% 
% 
% TransitionWidthDegDummy = CurvesIntersectionDeg - ThetaCusp - 1;
% TransitionWidthDeg = TransitionWidthDegDummy - mod(TransitionWidthDegDummy, DeltaThetaDeg); 
% TransitionWidth = TransitionWidthDeg/DeltaThetaDeg;

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
% rGuessBottom = InterpolatedGridInside(:);
rGuessBottom = InterpolatedGridInside(:);


%% Plotting the initial BOTTOM guess


figure;

    hold on
    % Ring
    Bp = 20000*10^(-9);     % Equatorial field, T
    Rp = 60280*10^3;        % Planet radius, m
    mu_0 = 4*pi*10^(-7);
    Pswnpa = 0.02;
    b0 = sqrt(2*mu_0*Pswnpa*10^(-9));   % Field scale
    r0 = (2*Bp*Rp^3/b0)^(1/3);          % Distance scale
    RpScaled = Rp/r0;
    r = 66900*10^3/r0; %inner radius
    R = 140180*10^3/r0; %outer radius
    steps = 100;
    theta = 0:2*pi/steps:2*pi;
    r_outer = ones(size(theta))*R;
    r_inner = ones(size(theta))*r;
    rr = [r_outer;r_inner];
    theta = [theta;theta];
    xx = rr.*cos(theta);
    yy = rr.*sin(theta);
    zz = zeros(size(xx));
    Rings = surf(xx, yy,zz);
    set(Rings,'FaceColor',[0.6000    0.5000    0.5000],...
              'FaceAlpha', 0.5,...
              'LineStyle','none');

    % CAN Disk
    a = 8*Rp/r0;      % Inner radius of the current distribution (in r0)
    b = 15.5*Rp/r0;   % Outer radius of the current distribution (in r0)
    D = 3*Rp/r0;       % Disk Semi-thickness (in r0)
    t = linspace(0,2*pi);
    rin = a;
    rout = b;
    center = [0, 0];
    xin = rin*cos(t);
    xout = rout*cos(t);
    yin = rin*sin(t);
    yout = rout*sin(t);
    z1 = -D;
    z2 = D;
    
    bottom = patch(center(1)+[xout,xin], ...
               center(2)+[yout,yin], ...
               z1*ones(1,2*length(xout)),'', 'FaceAlpha', 0.5, 'LineStyle', 'none');
    top = patch(center(1)+[xout,xin], ...
            center(2)+[yout,yin], ...
            z2*ones(1,2*length(xout)),'', 'FaceAlpha', 0.5, 'LineStyle', 'none');
    [X,Y,Z] = cylinder(1,length(xin));
    outer = surf(rout*X+center(1), ...
             rout*Y+center(2), ...
             Z*(z2-z1)+z1, 'FaceAlpha', 0.4, 'LineStyle', 'none');
    inner = surf(rin*X+center(1), ...
             rin*Y+center(2), ...
             Z*(z2-z1)+z1, 'FaceAlpha', 0.4, 'LineStyle', 'none');


    ThetaSpanBottom = Thetaq_Temp;
    PhiSpanBottom = Phiq_Temp;
    ConstructionSurfaceContours(InterpolatedGridTot.', ThetaSpanBottom, PhiSpanBottom);
    set(gcf,'color','w');
    axis([-1.5, 1.5, -2.3, 2.3, -2.3, 2.3])
    box on 
    view([45+90+5 30])
    grid on
    

%% Analysing Initial Bottom Grid

% theta = ThetaCusp+DeltaTheta : DeltaTheta : ThetaMax
% phi = DeltaPhi : DeltaPhi : PhiMax - DeltaPhi
NbPointsGridBottomInside = size(rGuessBottom,1) ;
NbPointsPhiBottomInside = size(BottomGrid, 1);
NbPointsThetaBottomInside = size(BottomGrid, 2);

BottomGridInitial = InterpolatedGridTot;
BottomGridInitialInside = BottomGridInitial(2:end-1, 2:end-1);
rTotalSurface = BottomGridInitialInside(:);
rBottom = BottomGridInitial(1,:).';
rTop = BottomGridInitial(end,:);
rLeft = BottomGridInitial(:,1);
figure;
PB_Assessment = F_values_Bottom(rTotalSurface, rBottom, rTop, rLeft, ThetaCusp, ThetaMaxDeg, PhiMaxDeg, DeltaThetaDeg, DeltaPhiDeg, NbPointsGridBottomInside, SystemParameters);
ErrorVector = log10(abs(PB_Assessment));
ErrorGrid = reshape(ErrorVector, length(rLeft)-2, length(rBottom)-2);
pcolor(ErrorGrid);
colorbar
caxis([-4 0])
    



%% Correcting the BOTTOM grid


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
rRight = ArrayGridTot(:,end);

rNew = fsolve( @(r) Values_Jacobian_Bottom( r, rBottom, rLeft, rTop, ThetaCusp, ThetaMaxDeg, PhiMaxDeg, DeltaThetaDeg, DeltaPhiDeg, NbPointsGridBottomInside, SystemParameters, epsilon), rGuessBottom, opts);

rTailRegionCorrected = rNew;


%% TEST


J = Jacobian_Numerical_Bottom(ThetaCusp, ThetaMaxDeg, PhiMaxDeg, DeltaThetaDeg, DeltaPhiDeg, rTailRegionCorrected, rLeft, rBottom, rTop, SystemParameters, epsilon);
pcolor(abs(J)); colorbar;

%% Plotting the corrected BOTTOM guess


figure;

    hold on
    % Ring
    Bp = 20000*10^(-9);     % Equatorial field, T
    Rp = 60280*10^3;        % Planet radius, m
    mu_0 = 4*pi*10^(-7);
    Pswnpa = 0.02;
    b0 = sqrt(2*mu_0*Pswnpa*10^(-9));   % Field scale
    r0 = (2*Bp*Rp^3/b0)^(1/3);          % Distance scale
    RpScaled = Rp/r0;
    r = 66900*10^3/r0; %inner radius
    R = 140180*10^3/r0; %outer radius
    steps = 100;
    theta = 0:2*pi/steps:2*pi;
    r_outer = ones(size(theta))*R;
    r_inner = ones(size(theta))*r;
    rr = [r_outer;r_inner];
    theta = [theta;theta];
    xx = rr.*cos(theta);
    yy = rr.*sin(theta);
    zz = zeros(size(xx));
    Rings = surf(xx, yy,zz);
    set(Rings,'FaceColor',[0.6000    0.5000    0.5000],...
              'FaceAlpha', 0.5,...
              'LineStyle','none');

    % CAN Disk
    a = 8*Rp/r0;      % Inner radius of the current distribution (in r0)
    b = 15.5*Rp/r0;   % Outer radius of the current distribution (in r0)
    D = 3*Rp/r0;       % Disk Semi-thickness (in r0)
    t = linspace(0,2*pi);
    rin = a;
    rout = b;
    center = [0, 0];
    xin = rin*cos(t);
    xout = rout*cos(t);
    yin = rin*sin(t);
    yout = rout*sin(t);
    z1 = -D;
    z2 = D;
    
    bottom = patch(center(1)+[xout,xin], ...
               center(2)+[yout,yin], ...
               z1*ones(1,2*length(xout)),'', 'FaceAlpha', 0.5, 'LineStyle', 'none');
    top = patch(center(1)+[xout,xin], ...
            center(2)+[yout,yin], ...
            z2*ones(1,2*length(xout)),'', 'FaceAlpha', 0.5, 'LineStyle', 'none');
    [X,Y,Z] = cylinder(1,length(xin));
    outer = surf(rout*X+center(1), ...
             rout*Y+center(2), ...
             Z*(z2-z1)+z1, 'FaceAlpha', 0.4, 'LineStyle', 'none');
    inner = surf(rin*X+center(1), ...
             rin*Y+center(2), ...
             Z*(z2-z1)+z1, 'FaceAlpha', 0.4, 'LineStyle', 'none');


    ThetaSpanBottom = Thetaq_Temp;
    PhiSpanBottom = Phiq_Temp;
    GridBottomCorrected = reshape(rTailRegionCorrected, length(Phiq_Temp)-2, length(Thetaq_Temp)-2);
    ConstructionSurfaceContours(GridBottomCorrected.', Thetaq_Temp(2:end-1), Phiq_Temp(2:end-1));
    set(gcf,'color','w');
    axis([-1.5, 1.5, -2.3, 2.3, -2.3, 2.3])
    box on 
    view([45+90+5 30])
    grid on
    




%% Merging the two sub-surfaces

BottomGridCorrected_temp1 = reshape(rNew, [NbPointsPhiBottomInside, NbPointsThetaBottomInside]);
BottomGridCorrected_temp2 = vertcat(rBottom(2:end-1).', BottomGridCorrected_temp1, rTop(2:end-1));
BottomGridCorrected_temp3 = horzcat(rLeft, BottomGridCorrected_temp2);

rRightCorrected = 2*BottomGridCorrected_temp3(:,end) - BottomGridCorrected_temp3(:,end-1);

BottomGridCorrected = horzcat(BottomGridCorrected_temp3, rRightCorrected);
% theta = rLeftBoundary : delta : ThetaMax (for delta = 1deg, 70:125)
% phi = 0:90

SurfaceCorrected = horzcat(SubSolarCapLong, BottomGridCorrected(:,2:end));

    
    
    ConstructionSurfaceContours(SurfaceCorrected.', ThetaSpan, PhiSpan);
    set(gcf,'color','w');
    axis([-1.5, 1.5, -2.3, 2.3, -2.3, 2.3])
    box on 
    view([45+90+5 30])
    grid on
%     
% %     
% %% Test: Merging with linear extrapolation
% 
% NighSidePosition = (90-ThetaCusp)/DeltaThetaDeg;
% BottomGridNightSide = SurfaceCorrected(:, NighSidePosition+1:end);
% BottomGridDaySide = SurfaceCorrected(:, 1:NighSidePosition);
% 
% NightSideTopBoundary = 2*BottomGridNightSide(end-1,:)-BottomGridNightSide(end-2,:);
% NightSideBottomBoundary = 2*BottomGridNightSide(3,:)-BottomGridNightSide(2,:);
% BottomGridNightSideExtrapolated = vertcat(NightSideBottomBoundary, BottomGridNightSide(2:end-1, :), NightSideTopBoundary);
% 
% SurfaceCorrected_Extrapolated = horzcat(BottomGridDaySide, BottomGridNightSideExtrapolated);
% 
% 
%     
% %     SurfaceCorrected = horzcat(SubSolarCapLong, BottomGridCorrected(:,1:end));
% 
%     
%     figure;
%     ConstructionSurfaceContours(SurfaceCorrected_Extrapolated.', ThetaSpan, PhiSpan);
%     set(gcf,'color','w');
%     axis([-1.5, 1.5, -2.3, 2.3, -2.3, 2.3])
%     box on 
%     view([45+90+5 30])
%     grid on
    
    
    
    
    
%% Analysing Corrected Bottom Grid

BottomGridCorrected = SurfaceCorrected(:,size(SubSolarCapLong, 2):end);
BottomGridInside = BottomGridCorrected(2:end-1, 2:end-1);
rTotalSurface = BottomGridInside(:);
rBottom = BottomGridCorrected(1,:).';
rTop = BottomGridCorrected(end,:);
rLeft = BottomGridCorrected(:,1);
figure;
PB_Assessment = F_values_Bottom(rTotalSurface, rBottom, rTop, rLeft, ThetaCusp, ThetaMaxDeg, PhiMaxDeg, DeltaThetaDeg, DeltaPhiDeg, NbPointsGridBottomInside, SystemParameters);
ErrorVector = log10(abs(PB_Assessment));
ErrorGrid = reshape(ErrorVector, length(rLeft)-2, length(rBottom)-2);
pcolor(ErrorGrid);
colorbar
caxis([-4 0])
    

%% Analysing the TOTAL surface

% SurfaceCorrectedTranspose = SurfaceCorrected.';
SurfaceCorrectedInside = SurfaceCorrected(2:end-1, 2:end-1);
rTotalSurface = SurfaceCorrectedInside(:);
rBottom = SurfaceCorrected(1,:);
rTop = SurfaceCorrected(end,:);
figure;
PB_Assessment = F_values_Numerical(rTotalSurface, rBottom, rTop, rSubSolarNose, ThetaMaxDeg, PhiMaxDeg, DeltaThetaDeg, DeltaPhiDeg, NbPointsGridInside, SystemParameters);
ErrorVector = log10(abs(PB_Assessment));
ErrorGrid = reshape(ErrorVector, length(PhiSpan)-2, length(ThetaSpan)-2);
pcolor(ErrorGrid);
colorbar
caxis([-4 0])



% %% plotting new equatorial and NMMP solutions
% figure;
% hold on
% 
% [Xequator, Yequator] = pol2cart(Thetaq_Temp.', rBottom);
% plot(Yequator, Xequator)
% 
% [Xmeridian, Ymeridian,] = pol2cart(Thetaq_Temp, rTop);
% plot(Ymeridian, Xmeridian)
% 
% [XequatorCurves, YequatorCurves] = pol2cart(ThetaSpanCurves.', rEquatorCurve);
% plot(YequatorCurves, XequatorCurves)
% 
% [XmeridianCurves, YmeridianCurves] = pol2cart(ThetaSpanCurves.', rMeridianCurve);
% plot(YmeridianCurves, XmeridianCurves)
% 
% axis equal
% axis([0 2.5 -1 1.5])
% grid on
% 
% hold off









%% Correction: Initialisation
epsilon = 1.e-8;

ArrayGridInside = SurfaceInitialTot(2:end-1, 2:end-1);
ArrayGridTot = SurfaceInitialTot;

rGuess = ArrayGridInside(:);

FValues = F_values_Numerical(rGuess, rEquator, rMeridian, rSubSolarNose, ThetaMaxDeg, PhiMaxDeg, DeltaThetaDeg, DeltaPhiDeg, NbPointsGridInside, SystemParameters);
FValuesGrid =  reshape(FValues, NbPointsThetaInside, NbPointsPhiInside);
figure;
pcolor(abs(FValuesGrid));
colorbar

%% Correction: Loop
NumberSteps = 5;

for k=1:NumberSteps

Jacobian =  NumericalJacobian(ThetaMaxDeg, PhiMaxDeg, DeltaThetaDeg, DeltaPhiDeg, rGuess, rSubSolarNose, rEquator, rMeridian, SystemParameters, epsilon);
rCorrection = -Jacobian\FValues;
rNew = rGuess + rCorrection;

FValuesCorrected = F_values_Numerical(rNew, rEquator, rMeridian, rSubSolarNose, ThetaMaxDeg, PhiMaxDeg, DeltaThetaDeg, DeltaPhiDeg, NbPointsGridInside, SystemParameters);
FValuesGridCorrected =  reshape(FValuesCorrected, NbPointsThetaInside, NbPointsPhiInside);
figure;
pcolor(abs(FValuesGridCorrected));
colorbar


rGuess = rNew;
FValues = FValuesCorrected;

end
pcolor(abs(Jacobian))
colorbar

ArrayCorrectedInside =  reshape(rNew, NbPointsPhiInside, NbPointsThetaInside);
ArrayCorrected = vertcat(rEquator(2:end-1).', ArrayCorrectedInside, rMeridian(2:end-1).');
ThetaSpanInside = ThetaSpan(2:end-1);
PhiSpanInside = PhiSpan(1:end);
figure;
ConstructionSurface(ArrayCorrected.',ThetaSpanInside,PhiSpanInside)


%% Correction: fsolve with Jacobian


NumIterations = 5;

opts = optimoptions('fsolve','algorithm','Levenberg-Marquardt', 'Display', 'iter');
opts.MaxIterations = NumIterations;
opts.MaxFunEvals = 10;
opts.SpecifyObjectiveGradient = true;
opts.CheckGradients = false;
% opts.JacobPattern = JacobianPattern(ThetaMaxDegSubSolar, PhiMaxDegSubSolar, DeltaThetaDeg, DeltaPhiDeg);
opts.FiniteDifferenceType = 'central';
rss = rSubSolarNose;
[rNew] = fsolve( @(r) F_valuesJacobian_Numerical(r, rEquator, rMeridian, rSubSolarNose, ThetaMaxDeg, PhiMaxDeg, DeltaThetaDeg, DeltaPhiDeg, NbPointsGridInside, SystemParameters, epsilon), rGuess, opts);


ArrayCorrectedInside =  reshape(rNew, NbPointsPhiSubSolarInside, NbPointsThetaSubSolarInside);
ArrayCorrected = vertcat(rEquator(2:end-1).', ArrayCorrectedInside, rMeridian(2:end-1).');
ThetaSpanInside = (DeltaThetaDeg:DeltaThetaDeg:ThetaMaxDegSubSolar-DeltaThetaDeg)*pi/180;
PhiSpanInside = (DeltaPhiDeg:DeltaPhiDeg:PhiMaxDegSubSolar-DeltaThetaDeg)*pi/180;

figure;
ConstructionSurface(ArrayCorrectedInside.',ThetaSpanInside,PhiSpanInside)
%}
