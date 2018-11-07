function [SurfaceInitialTot, ThetaCusp, rEquator, rMeridian] = InitialGuess(ThetaMaxDeg, PhiMaxDeg, DeltaThetaDeg, DeltaPhiDeg, Plots)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs
%    - ThetaMaxDeg, PhiMaxDeg: scalars, maximum values for theta and phi,
%           in degrees
%    - DeltaThetaDeg, DeltaPhiDeg: scalars, angular increments, in degrees
%    - Plots: list of strings, selection of plots to display
%           - "CurveEquator": solution in the equatorial plane
%           - "CurveMeridianParts: solution in the noon-midnight meridian
%               plane, 2 curves
%           - "CurveMeridian": solution in the noon-midnight meridian,
%               merged
%           - "SurfaceInitialTot": initial 3D surface
%           - "GridAnalysis": colourmap assessing the pressure balance on
%               the planar grid
%           - "GridWrapped": colourmap assessing the pressure balance
%           projected onto the 3D surface
%
% Outputs
%     - SurfaceInitialTot: phi*theta array, initial surface
%     - ThetaCusp: scalar, theta value of the cusp, in degrees
%     - rEquator: vector, solution in the equatorial plane
%           theta = 0:DeltaThetaDeg:ThetaMaxDeg
%     - rMeridian: vector, solution in the noon-midnight meridian plane
%           theta = 0:DeltaThetaDeg:ThetaMaxDeg
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Definition of the grid

ThetaSpan = (0 : DeltaThetaDeg : ThetaMaxDeg)*pi/180;
PhiSpan = (0 : DeltaPhiDeg: PhiMaxDeg).*pi/180;

NbPointsTheta = ThetaMaxDeg/DeltaThetaDeg+1; % Number of points in theta-direction for phi=cst (theta = 0:Max_deg_equ)
NbPointsPhi = PhiMaxDeg/DeltaPhiDeg+1; % Number of points in phi-direction for theta=cst (phi = 0:Max_deg_phi)
NbPointsGrid = NbPointsTheta*NbPointsPhi;

NbPointsThetaInside = NbPointsTheta-2; % Number of points on the equator without first (theta=0) and last (theta=theta_max) points
NbPointsPhiInside = NbPointsPhi-2;
NbPointsGridInside = NbPointsThetaInside*NbPointsPhiInside; % Number of points on the interior grid for theta, phi = delta : Max_deg-delta



%% Equatorial plane: explicit integration

[ThetaEquator, rEquator] = Balance_Equator(ThetaSpan);
% theta = 0 : delta_theta : Max_deg_equ
% Nb_points_theta = Max_deg_equ/delta_theta_deg+1;



%% Meridian plane: piece-wise integration

%%% Meridian Plane: solving explicitely for the PRE-cusp section
r0 = 1; % Standoff distance / distance scale
rSubSolarNose = 1;
ThetaSpanMeridian = [0.001, DeltaThetaDeg:DeltaThetaDeg:ThetaMaxDeg]*pi/180;
funPreCusp = @(theta, r) drMeridian(theta, r);

[ThetaMeridian1, rMeridian1] = ode45( funPreCusp, ThetaSpanMeridian, r0 );


%%% Meridian Plane: solving explicitely for the POST-cusp section, in 2 directions
r0 = 2^(1/3); % Point A2 at theta = pi/2 in meridian plane
thetaSpanAfter = (90:DeltaThetaDeg:ThetaMaxDeg)*(pi/180);
funPostCusp = @(theta, r) drMeridian2(theta, r);
[thetaAfter, rAfter] = ode45( funPostCusp, thetaSpanAfter, r0 );

thetaSpanPrevious = (90:-DeltaThetaDeg:DeltaThetaDeg)*pi/180;
[thetaPrevious, rPreviousTemp] = ode45( funPostCusp, thetaSpanPrevious, r0 );

rPrevious = [0; flipud(rPreviousTemp(2:end))];
% Segment upstream of A2


%%% Meridian Plane: Merging the two curves
ThetaMerged = horzcat( 0, flipud(thetaPrevious(1:1:end)).', thetaSpanAfter(2:end));
rTwoCurves = [rPrevious; rAfter];

rMerged = rTwoCurves.';

C=0;
for k = 1:length(ThetaMerged)
    
    if rMerged(k) < rMeridian1(k)
       rMerged(k) = rMeridian1(k);
       C = C+1;
    end
    
end

ThetaCusp = DeltaThetaDeg * ( C-1 ) ;

rMeridian = rMerged;



%% Construction of the initial surface

Array = zeros(size(ThetaSpan,2),size(PhiSpan,2));       % Initialising the Array
Array(1,:) = ones(size(PhiSpan));                       % Left-boundary = sub-solar nose
Array(:,1) = rEquator;                                  % Bottom boundary = equator
Array(:,end) = rMeridian;                               % Upper boundary = meridian


PhiSpank =  (0 : DeltaPhiDeg : PhiMaxDeg).*pi/180 ;

for k = 2:size(ThetaSpan,2)
    for m = 1:size(PhiSpank,2)
        a = rEquator(k)*sin(ThetaSpan(k)); 
        b = rMeridian(1*k)*sin(ThetaSpan(k));
        c = sqrt(a^2-b^2);
        phi = PhiSpank(m);
        zk = (a*b) / sqrt(b^2*cos(phi)^2+a^2*sin(phi)^2);
        Array(k,m) = ( (a*b) / sqrt(b^2*cos(phi)^2+a^2*sin(phi)^2) ) / sin(ThetaSpan(k)) ;
    end
end

% Correction of the 'theta_max' row, via linear extrapolation

for k = 1:size(PhiSpan,2)

    Array(end, k) = 2*Array(end-1, k) - Array(end-2, k);
    
end

SurfaceInitialTot = Array.';



%% Analysis of the initial surface

ArrayGridInside = Array(2:end-1, 2:end-1).';
rGuess = ArrayGridInside(:);





%% Plots
% Curve solution in the equatorial plane
if ismember("CurveEquator", Plots)
    
    figure;
    title('Solution in the equatorial plane')
    hold on
        ThetaAxis = pi/2-ThetaEquator;
        [Yeq, Zeq] = pol2cart(ThetaAxis, rEquator);
        EquatorialSolution = plot(Yeq, Zeq, 'LineWidth', 4);    
        set(EquatorialSolution, 'Color', [00.8500    0.6000    0.4000]); 
        
        axes = gca;
        set(axes, 'FontName'   , 'Helvetica',...
                  'XAxisLocation', 'origin',...
                  'YAxisLocation', 'origin',...
                  'FontSize', 16,...
                  'LineWidth'   , 1 ...
              );
        Ylabel = ylabel('Z (r_0)');
        Xlabel = xlabel('X (r_0)');  
        set([Xlabel, Ylabel], 'FontName'   , 'AvantGarde', 'FontSize', 18);
        set(Ylabel, 'position', [-0.42 1.05 0]);
        axis equal
        axis([-0.5 1.7 -0.9 1.1]);
        grid on
        box on
        set(gcf,'color','w')
        
        % Annotation: Solar wind and Magnetic moment
        c_SolarWind = [0.6400    0.0800    0.1800];
        c_MagneticMoment = [0.6400    0.0800    0.1800];
        mArrow3([1.3,0.90,0],[1.3,0.7,0],'color',c_SolarWind,'stemWidth', 0.0085, 'tipWidth', 0.018, 'FaceAlpha', 0.9);
        SolarWind = '$\vec{v}$';
        MagneticMomentVector = '$\vec{\mathcal{M}}$';
        text(1.35, 0.8, SolarWind, 'Interpreter','latex', 'FontSize', 20, 'FontName', 'Helvetica', 'color', c_SolarWind);
        rectangle('Position',[-0.06 -0.06 0.12 0.12],'Curvature',[1 1], 'FaceColor', c_MagneticMoment, 'Edgecolor', 'none')
        rectangle('Position',[-0.05 -0.05 0.05*2 0.05*2],'Curvature',[1 1], 'FaceColor', 'white', 'Edgecolor', 'none')
        text(0.1, -0.1, MagneticMomentVector, 'Interpreter','latex', 'FontSize', 16, 'FontWeight', 'bold', 'FontName', 'Helvetica', 'color', c_MagneticMoment);
        line([-0.05, 0.05], [0, 0], 'LineWidth', 1, 'Color', [0.1500    0.1500    0.1500]);
        line([0, 0], [-0.05, 0.05], 'LineWidth', 1, 'Color', [0.1500    0.1500    0.1500]);
        rectangle('Position',[-0.017 -0.017 0.017*2 0.017*2],'Curvature',[1 1], 'FaceColor', c_MagneticMoment, 'Edgecolor', 'none');
        
        % Annotation: Coordinate system
        c_Position = [0.3000    0.4500    0.7400];
        line([0,Yeq(100)], [0, Zeq(100)], 'LineWidth', 2, 'LineStyle', ':', 'Color', c_Position);
        circular_arrow(EquatorialSolution, 0.3, [0,0], 65, 30, 1, c_Position, 7, 'vback2', 2)
        PositionTheta = '$\theta$';
        PositionR = '$r$';
        text(0.15, 0.34, PositionTheta, 'Interpreter','latex', 'FontSize', 20, 'FontWeight', 'bold', 'FontName', 'Helvetica', 'color', c_Position);
        text(0.53, 0.35, PositionR, 'Interpreter','latex', 'FontSize', 22, 'FontWeight', 'bold', 'FontName', 'Helvetica', 'color', c_Position);
	
	hold off
    
end



%% 
% Piecewise curve solutions in the noon-midnight meridian plane

if ismember("CurveMeridianParts", Plots)
    
    figure;
    title('Piece-wise solution in the noon-midnight meridian plane')
    hold on
        [Y, Z] = pol2cart(pi/2-ThetaMeridian1, rMeridian1);
        MeridianCurve1 = plot(Y, Z, 'LineWidth', 4);
        [Y, Z] = pol2cart(pi/2-ThetaMerged, rTwoCurves.');
        MeridianCurve2 = plot(Y, Z, 'LineWidth', 4);
        
        set(MeridianCurve1, 'Color', [0.4000    0.4500    0.7400]); 
        set(MeridianCurve2, 'Color', [0.4700    0.6700    0.5000]); 
        
        axes = gca;
        set(axes, 'FontName'   , 'Helvetica',...
                  'XAxisLocation', 'origin',...
                  'YAxisLocation', 'origin',...
                  'FontSize', 16,...
                  'LineWidth'   , 1 ...
              );
        Ylabel = ylabel('Z (r_0)');
        Xlabel = xlabel('Y (r_0)');  
        set([Xlabel, Ylabel], 'FontName'   , 'AvantGarde', 'FontSize', 18);
        set(Ylabel, 'position', [-0.42 1.05 0]);
        axis equal
        axis([-0.5 1.7 -0.9 1.1]);
        grid on
        box on
        set(gcf,'color','w')        
        
  
        % Annotation: Solar wind and Magnetic moment
        c_SolarWind = [0.6400    0.0800    0.1800];
        c_MagneticMoment = [0.6400    0.0800    0.1800];
        mArrow3([1.3,0.90,0],[1.3,0.7,0],'color',c_SolarWind,'stemWidth', 0.0085, 'tipWidth', 0.018, 'FaceAlpha', 0.9);
        SolarWind = '$\vec{v}$';
        MagneticMomentVector = '$\vec{\mathcal{M}}$';
        text(1.35, 0.8, SolarWind, 'Interpreter','latex', 'FontSize', 20, 'FontName', 'Helvetica', 'color', c_SolarWind);

        mArrow3([0.18, 0, 0], [-0.18, 0, 0],'color',c_MagneticMoment,'stemWidth', 0.01, 'tipWidth', 0.02);  
        text(0.1, -0.1, MagneticMomentVector, 'Interpreter','latex', 'FontSize', 16, 'FontWeight', 'bold', 'FontName', 'Helvetica', 'color', c_MagneticMoment);

        % Annotation: Coordinate system
        c_Points = [0.7000    0.3300    0.1000];
        A1 = plot(1, 0, 'o', 'color', c_Points, 'MarkerSize', 10, 'LineWidth', 2);
        A2 = plot(2^(1/3), 0, 'o', 'color', c_Points, 'MarkerSize', 10, 'LineWidth', 2);
        text(1-0.17, -0.07, '$A_1$', 'Interpreter','latex', 'FontSize', 20, 'FontName', 'Helvetica', 'color', c_Points);
        text(2^(1/3)-0.15, -0.07, '$A_2$', 'Interpreter','latex', 'FontSize', 20, 'FontName', 'Helvetica', 'color', c_Points);
        
  	hold off
    
end



%%
% Curve solution in the noon-midnight meridian plane

if ismember("CurveMeridian", Plots)
    
    figure;
    title('Solution in the noon-midnight meridian plane')
    hold on
        [Y, Z] = pol2cart(pi/2-ThetaMerged, rMerged);
        MeridianSolution = plot(Y, Z, 'LineWidth', 4);
        
        PointXCusp = [0, rSubSolarNose*sin( ThetaCusp*pi/180 )];
        PointYCusp = [0  rSubSolarNose*cos( ThetaCusp*pi/180 )]; 
        c_Position = [0.4900    0.1800    0.5000];
        CuspPosition = plot(PointXCusp, PointYCusp, ':', 'LineWidth', 2, 'Color', c_Position);
        C = plot(PointXCusp(2), PointYCusp(2), 'o', 'color', c_Position, 'MarkerSize', 10, 'LineWidth', 2);

        set(MeridianSolution, 'Color', [0.4000    0.6000    0.3000]); 
        
        axes = gca;
        set(axes, 'FontName'   , 'Helvetica',...
                  'XAxisLocation', 'origin',...
                  'YAxisLocation', 'origin',...
                  'FontSize', 16,...
                  'LineWidth'   , 1 ...
              );
        Ylabel = ylabel('Z (r_0)');
        Xlabel = xlabel('Y (r_0)');  
        set([Xlabel, Ylabel], 'FontName'   , 'AvantGarde', 'FontSize', 18);
        set(Ylabel, 'position', [-0.42 1.05 0]);
        axis equal
        axis([-0.5 1.7 -0.9 1.1]);
        grid on
        box on
        set(gcf,'color','w')
        

        % Annotation: Solar wind and Magnetic moment
        c_SolarWind = [0.6400    0.0800    0.1800];
        c_MagneticMoment = [0.6400    0.0800    0.1800];
        mArrow3([1.3,0.90,0],[1.3,0.7,0],'color',c_SolarWind,'stemWidth', 0.0085, 'tipWidth', 0.018, 'FaceAlpha', 0.9);
        SolarWind = '$\vec{v}$';
        MagneticMomentVector = '$\vec{\mathcal{M}}$';
        text(1.35, 0.8, SolarWind, 'Interpreter','latex', 'FontSize', 20, 'FontName', 'Helvetica', 'color', c_SolarWind);

        mArrow3([0.18, 0, 0], [-0.18, 0, 0],'color',c_MagneticMoment,'stemWidth', 0.01, 'tipWidth', 0.02);  
        text(0.1, -0.1, MagneticMomentVector, 'Interpreter','latex', 'FontSize', 16, 'FontWeight', 'bold', 'FontName', 'Helvetica', 'color', c_MagneticMoment);
        
        % Annotation: Coordinate system
        c_Points = [0.7000    0.3300    0.1000];
        circular_arrow(MeridianSolution, 0.3, [0,0], 55, 50, 1, c_Position, 7, 'vback2', 2);
        PositionTheta = '$\theta_C$';
        text(0.2, 0.32, PositionTheta, 'Interpreter','latex', 'FontSize', 20, 'FontWeight', 'bold', 'FontName', 'Helvetica', 'color', c_Position);
        A2 = plot(2^(1/3), 0, 'o', 'color', c_Points, 'MarkerSize', 10, 'LineWidth', 2);
        text(2^(1/3)-0.15, -0.07, '$A_2$', 'Interpreter','latex', 'FontSize', 20, 'FontName', 'Helvetica', 'color', c_Points);
        text(PointXCusp(2)+0.05, PointYCusp(2)+0.05, '$C$', 'Interpreter','latex', 'FontSize', 20, 'FontName', 'Helvetica', 'color', c_Position);

	hold off
    
end



%%
% Initial surface to be used as an initial 'guess'

if ismember('SurfaceInitialTot', Plots)
    
    figure;
    title('Initial Surface')
    
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
    surf(yy,zz,xx, 'FaceAlpha', 0.3)
    
    ConstructionSurface(Array,ThetaSpan,PhiSpan)
    set(gcf,'color','w');
    axis([-1.7, 1.7, -1.7, 1.7, -1.7, 1])
    box on 
    grid on
    
end



%%
% Iso-theta contours of the initial surface

if ismember('SurfaceInitialTotContours', Plots)
    
    figure;
    title('Initial Surface contours')
    
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
    surf(yy,zz,xx, 'FaceAlpha', 0.3)
    
    ConstructionSurfaceContours(Array,ThetaSpan,PhiSpan)
    set(gcf,'color','w');
    axis([-1.7, 1.7, -1.7, 1.7, -1.7, 1])
    box on 
    grid on
    
end



%%
% Initial surface to be used as an initial 'guess' + Field lines

if ismember("SurfaceInitialTotFieldLines", Plots)

    figure;
    title('Initial surface and dipole field lines')  

        L = (0.1:0.5:50).';
    theta = (0:0.1:180)*pi/180;
    phi = 90*pi/180;

    hold on

    
    for m = 1:length(phi)
      for k = 1:length(L)
         r = L(k).*sin(theta).^2;
         [X, Y, Z] = sph2cart(phi(m), theta, r);
         axis([-1.7, 1.7, -1.7, 1.7, -1.7, 1.7]);
         plot3(X,Y,Z, ':', 'LineWidth', 2)
      end
    end

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
    surf(yy,zz,xx, 'FaceAlpha', 0.3)
    
    ConstructionSurface(Array,ThetaSpan,PhiSpan)
    set(gcf,'color','w');
    box on
    grid on
    
end



%%
% Absolute pressure balance assessment on the 2D grid

if ismember('GridAnalysis', Plots)
    
    figure;
    title('Absolute PB assessment on the the initial surface')  
    FValuesGridInsideVector = F_values(rGuess, rEquator, rMeridian, 1, ThetaMaxDeg, PhiMaxDeg, DeltaThetaDeg, DeltaPhiDeg, NbPointsGridInside);
    FValuesGridInside = reshape(FValuesGridInsideVector,[NbPointsPhiInside, NbPointsThetaInside]);
    ErrorAbsolute = (abs(FValuesGridInside(1:1:end,1:1:end)));
    pcolor(ErrorAbsolute(1:end, 1:end)),shading flat; colorbar; 
        title('Pressure balance assessment on the initial surface')
        axes = gca;
        xlabel('Theta','FontSize',16,'FontWeight','bold');
        xticks(0:10/DeltaThetaDeg:NbPointsThetaInside);
        xticklabels({(0:10/DeltaThetaDeg:NbPointsThetaInside)*DeltaThetaDeg});
        ylabel('Phi','FontSize',16,'FontWeight','bold');
        yticks(0:10/DeltaPhiDeg:NbPointsPhiInside);
        yticklabels({(0:10/DeltaPhiDeg:NbPointsPhiInside)*DeltaPhiDeg});
        set(gca, 'FontSize', 16);
        axes.XAxisLocation = 'origin';
        axes.YAxisLocation = 'origin';
        set(gcf,'color','w');
        box on
        
end



%%
% Absolute pressure balance assessment projected onto the 3D intial surface

if ismember('GridWrappedAbsolute', Plots)
    
    figure();
    FValuesGridInsideVector = F_values(rGuess, rEquator, rMeridian, 1, ThetaMaxDeg, PhiMaxDeg, DeltaThetaDeg, DeltaPhiDeg, NbPointsGridInside);
    FValuesGridInside = reshape(FValuesGridInsideVector,[NbPointsPhiInside, NbPointsThetaInside]);
    ErrorAbsolute = (abs(FValuesGridInside(1:1:end,1:1:end)));
    pcolor(ErrorAbsolute(1:end, 1:end));shading flat;
    MaxErrorScale = 0.15;
        caxis([0 MaxErrorScale]);
        axes = gca;
        xlabel([]);
        xticks([]);
        xticklabels([]);
        ylabel([]);
        yticks([]);
        yticklabels([]);
        set(gca, 'FontSize', 16);
        axes.XAxisLocation = 'origin';
        axes.YAxisLocation = 'origin';
        set(gcf,'color','w');
        box on
        camroll(-90)
        export_fig ColourMap.png
    close();

    figure;
    title('Absolute PB assessment projected onto the initial surface')  
    
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
    surf(yy,zz,xx, 'FaceAlpha', 0.3)
    
    ImgRGB = imread('ColourMap.png');
    [ImRows, ImCols, ImPlanes] = size(ImgRGB);
    ConstructionSurfaceWarp(Array(1:1:end, 1:1:end),ThetaSpan(1:1:end), PhiSpan(1:1:end), ImgRGB)
    set(gcf,'color','w');
%     colorbar('southoutside');
    colorbar('southoutside', 'YTickLabel', {'0.01', '', '0.1', '', '1', '', '10', '', '100'});
    caxis([0 MaxErrorScale]);
    view([0 90])
    
end



%%
% Absolute pressure balance assessment projected onto the 3D intial surface
% + Field lines

if ismember("GridWrappedAbsoluteFieldLines", Plots)

    L = (0.1:0.5:20).';
    theta = (0:0.1:180)*pi/180;
    phi = 90*pi/180;

    NbPointsTheta = ThetaMaxDeg/DeltaThetaDeg+1;   % Number of points in theta-direction for phi=cst (theta = 0:Max_deg_equ)
    NbPointsPhi = PhiMaxDeg/DeltaPhiDeg+1;         % Number of points in phi-direction for theta=cst (phi = 0:Max_deg_phi)
    NbPointsGrid = NbPointsTheta*NbPointsPhi;
    NbPointsThetaInside = NbPointsTheta-2;                           % Number of points on the equator without first (theta=0) and last (theta=theta_max) points
    NbPointsPhiInside = NbPointsPhi-2;
    NbPointsGridInside = NbPointsThetaInside*NbPointsPhiInside;      % Number of points on the interior grid for theta, phi = delta : Max_deg-delta

    PhiSpanConcatenated = (0:DeltaPhiDeg:PhiMaxDeg)*pi/180;
    ThetaSpanConcatenated = (0:DeltaThetaDeg:ThetaMaxDeg)*pi/180;

    % Texture
    figure();
    FValuesGridInsideVector = F_values(rGuess, rEquator, rMeridian, 1, ThetaMaxDeg, PhiMaxDeg, DeltaThetaDeg, DeltaPhiDeg, NbPointsGridInside);
    FValuesGridInside = reshape(FValuesGridInsideVector,[NbPointsPhiInside, NbPointsThetaInside]);
    ErrorAbsolute = (abs(FValuesGridInside(1:1:end,1:1:end)));
    pcolor(ErrorAbsolute(1:end, 1:end));shading flat;
    MaxErrorScale = 0.15;
        caxis([0 MaxErrorScale]);
        axes = gca;
        xlabel([]);
        xticks([]);
        xticklabels([]);
        ylabel([]);
        yticks([]);
        yticklabels([]);
        set(gca, 'FontSize', 16);
        axes.XAxisLocation = 'origin';
        axes.YAxisLocation = 'origin';
        set(gcf,'color','w');
        box on
        camroll(-90)
        export_fig ColourMap.png
    close();

    figure;
    hold on
    
    % L-shells
    for m = 1:length(phi)
      for k = 1:length(L)
         r = L(k).*sin(theta).^2;
         [X, Y, Z] = sph2cart(phi(m), theta, r);
         axis([-1.7, 1.7, -1.7, 1.7, -1.7, 1.7]);
         plot3(X,Y,Z, ':', 'LineWidth', 2)
      end
    end
    
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
    surf(yy,zz,xx, 'FaceAlpha', 0.3)
    
    % Surface
    title('Absolute PB assessment projected onto the initial surface')  
    ImgRGB = imread('ColourMap.png');
    [ImRows, ImCols, ImPlanes] = size(ImgRGB);
    ConstructionSurfaceWarp(Array(1:1:end, 1:1:end),ThetaSpan(1:1:end), PhiSpan(1:1:end), ImgRGB)
    set(gcf,'color','w');
%     colorbar('southoutside');
    colorbar('southoutside', 'YTickLabel', {'0.01', '', '0.1', '', '1', '', '10', '', '100'});
    caxis([0 MaxErrorScale]);
    view([0 90])
    axis([-1.7, 1.7, -1.7, 1.7, -1.7, 1.7])
    box on
    grid on

end



%%
% Relative pressure balance assessment projected onto the 3D intial surface

if ismember('GridWrappedRelative', Plots)
  
    figure();
    FValuesGridInsideVector = F_valuesRelative(rGuess, rEquator, rMeridian, 1, ThetaMaxDeg, PhiMaxDeg, DeltaThetaDeg, DeltaPhiDeg, NbPointsGridInside);
    FValuesGridInside = reshape(FValuesGridInsideVector,[NbPointsPhiInside, NbPointsThetaInside]);
    ErrorAbsolute = log10(abs(FValuesGridInside(1:1:end,1:1:end)));
    pcolor(ErrorAbsolute(1:end, 1:end));shading flat;
    MaxErrorScale = 0.50;
        caxis([-4 0]);
        axes = gca;
        xlabel([]);
        xticks([]);
        xticklabels([]);
        ylabel([]);
        yticks([]);
        yticklabels([]);
        set(gca, 'FontSize', 16);
        axes.XAxisLocation = 'origin';
        axes.YAxisLocation = 'origin';
        set(gcf,'color','w');
        box on
        camroll(-90)
        export_fig ColourMap.png
    close();

    figure;
    title('Relative PB assessment projected onto the initial surface')  
    
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
    hold on
    Rings = surf(xx, yy,zz);
    set(Rings,'FaceColor',[0.6000    0.5000    0.5000],...
              'FaceAlpha', 0.5,...
              'LineStyle','none');
hold off
    
    
    ImgRGB = imread('ColourMap.png');
    [ImRows, ImCols, ImPlanes] = size(ImgRGB);
    ConstructionSurfaceWarp(Array(1:1:end, 1:1:end),ThetaSpan(1:1:end), PhiSpan(1:1:end), ImgRGB)
    axis([-1.7 1.7 -1.7 1.7 -1.7 1.7])
    CB = colorbar('eastoutside', 'YTickLabel', {'0.01', '', '0.1', '', '1', '', '10', '', '100'}, 'FontSize', 16, 'FontName', 'Helvetica');    
    CBPosition = CB.Position;
    set(CB, 'LineWidth', 1, 'Position', [CBPosition(1)+0.035 CBPosition(2)+0.05 0.03 0.67])
    caxis([-4 0]);
    view([45+90+5 20])
        
end



%%
% Relative pressure balance assessment projected onto the 3D intial surface
% + Field lines

if ismember("GridWrappedRelativeFieldLines", Plots)

    L = (0.1:0.5:20).';
    theta = (0:0.1:180)*pi/180;
    phi = 90*pi/180;

    NbPointsTheta = ThetaMaxDeg/DeltaThetaDeg+1;   % Number of points in theta-direction for phi=cst (theta = 0:Max_deg_equ)
    NbPointsPhi = PhiMaxDeg/DeltaPhiDeg+1;         % Number of points in phi-direction for theta=cst (phi = 0:Max_deg_phi)
    NbPointsGrid = NbPointsTheta*NbPointsPhi;
    NbPointsThetaInside = NbPointsTheta-2;                           % Number of points on the equator without first (theta=0) and last (theta=theta_max) points
    NbPointsPhiInside = NbPointsPhi-2;
    NbPointsGridInside = NbPointsThetaInside*NbPointsPhiInside;      % Number of points on the interior grid for theta, phi = delta : Max_deg-delta

    PhiSpanConcatenated = (0:DeltaPhiDeg:PhiMaxDeg)*pi/180;
    ThetaSpanConcatenated = (0:DeltaThetaDeg:ThetaMaxDeg)*pi/180;

    % Texture
    figure();
    FValuesGridInsideVector = F_valuesRelative(rGuess, rEquator, rMeridian, 1, ThetaMaxDeg, PhiMaxDeg, DeltaThetaDeg, DeltaPhiDeg, NbPointsGridInside);
    FValuesGridInside = reshape(FValuesGridInsideVector,[NbPointsPhiInside, NbPointsThetaInside]);
    ErrorAbsolute = log10(abs(FValuesGridInside(1:1:end,1:1:end)));
    pcolor(ErrorAbsolute(1:end, 1:end));shading flat;
    MaxErrorScale = 0.50;
        caxis([-4 0]);
        axes = gca;
        xlabel([]);
        xticks([]);
        xticklabels([]);
        ylabel([]);
        yticks([]);
        yticklabels([]);
        set(gca, 'FontSize', 16);
        axes.XAxisLocation = 'origin';
        axes.YAxisLocation = 'origin';
        set(gcf,'color','w');
        box on
        camroll(-90)
        export_fig ColourMap.png
    close();

    figure;
    hold on
    
    % L-shells
    for m = 1:length(phi)
      for k = 1:length(L)
         r = L(k).*sin(theta).^2;
         [X, Y, Z] = sph2cart(phi(m), theta, r);
         axis([-1.7, 1.7, -1.7, 1.7, -1.7, 1.7]);
         plot3(Z, X,Y, ':', 'LineWidth', 2)
      end
    end
    
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
    surf(xx, yy,zz, 'FaceAlpha', 0.3)
    
    % Surface
    title('Relative PB assessment projected onto the initial surface')  
    ImgRGB = imread('ColourMap.png');
    [ImRows, ImCols, ImPlanes] = size(ImgRGB);
    ConstructionSurfaceWarp(Array(1:1:end, 1:1:end),ThetaSpan(1:1:end), PhiSpan(1:1:end), ImgRGB)
    axis([-1.7 1.7 -1.7 1.7 -1.7 1.7])
    CB = colorbar('eastoutside', 'YTickLabel', {'0.01', '', '0.1', '', '1', '', '10', '', '100'}, 'FontSize', 16, 'FontName', 'Helvetica');    
    CBPosition = CB.Position;
    set(CB, 'LineWidth', 1, 'Position', [CBPosition(1)+0.035 CBPosition(2)+0.05 0.03 0.67])
    caxis([-4 0]);
    view([45+90+5 20])

end


end