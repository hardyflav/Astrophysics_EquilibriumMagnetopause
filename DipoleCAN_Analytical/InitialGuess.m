function [SurfaceInitialTot, ThetaCusp, rEquator, rMeridian, rSubSolarNose] = InitialGuess(SystemParameters, ThetaMaxDeg, PhiMaxDeg, DeltaThetaDeg, DeltaPhiDeg, Plots)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs
%    - SystemParameters: structure with parameters related to the planet,
%            hot plasma pressure and pressure balance
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


%% System Definition
beta = SystemParameters.beta;
M = SystemParameters.M;

%% Grid Definition

[~, NbPointsGridInside,...
NbPointsThetaInside, NbPointsPhiInside,...
~, ~] = GridDetails(ThetaMaxDeg, PhiMaxDeg, DeltaThetaDeg, DeltaPhiDeg);
ThetaSpan = (0 : DeltaThetaDeg : ThetaMaxDeg)*pi/180;
PhiSpan = (0 : DeltaPhiDeg: PhiMaxDeg).*pi/180;



%% Equatorial plane: explicit integration

% Finding the position of the sub-solar nose
thetaSubSolarNose = 0;
phiSubSolarNose = 0;
drdthetaSubSolarNose = 0;
drdphiSubSolarNose = 0;

% CAN disk parameters
CAN_DiskParameters = SystemParameters.CAN_DiskParameters;

% Solving the PB on the (OZ)-axis to find the subsolar nose
PB = @(r) PressureBalanceCANDisk(thetaSubSolarNose, phiSubSolarNose, r, drdthetaSubSolarNose, drdphiSubSolarNose, beta, SystemParameters);
rGuess = 1;
rSubSolarNose = fzero(PB, rGuess);


% Solving the PB in the equatorial plane from the subsolar nose
r_ini = rSubSolarNose;                                        % Initial "nose" distance, r0
dr_ini = 0;
opts = odeset('RelTol',1e-12,'AbsTol',1e-12);
ThetaSpanEquator = (0:DeltaThetaDeg:ThetaMaxDeg)*pi/180;
[ThetaEquator, rEquator] =  ode15i( @(theta, r, dr) PressureBalanceCANDisk(theta, 0, r, dr, 0, beta, SystemParameters), ThetaSpanEquator, r_ini, dr_ini, opts );
% theta = 0 : delta_theta : Max_deg_equ
% Nb_points_theta = Max_deg_equ/delta_theta_deg+1;
%  
% figure;
% hold on
% [X,Y] = pol2cart(ThetaEquator, rEquator);
% plot(Y, X, 'LineWidth', 3)
% grid on
% axis equal
% % axis([0, 1.5, -0.9, 1])
%         axes = gca;
%         set(gca, 'FontSize', 16)
%         axes.XAxisLocation = 'origin';
%         axes.YAxisLocation = 'origin';

%         
% %% TEST - Solving in the Night Side separately
% ThetaXaxis = pi/2;
% PhiXaxis = 0;
% rXaxis = rEquator(90/DeltaThetaDeg +1);
% drdphiXaxis = 0;
% dr_ini_Before = (rEquator(90/DeltaThetaDeg +1) - rEquator(90/DeltaThetaDeg +1-1))/(ThetaEquator(90/DeltaThetaDeg +1)-ThetaEquator(90/DeltaThetaDeg +1-1));
% 
% PB_After = @(drdtheta) PressureBalanceCANDisk(ThetaXaxis, PhiXaxis, rXaxis, drdtheta, drdphiXaxis, beta, SystemParameters);
% drdthetaXaxis = fsolve(PB_After, dr_ini_Before);
% 
% r_ini_After = rXaxis;
% dr_ini_After = drdthetaXaxis;
% ThetaSpanEquatorAfter = (90:DeltaThetaDeg:120)*pi/180;
% [ThetaEquatorAfter, rEquatorAfter] =  ode15i( @(theta, r, dr) PressureBalanceCANDisk(theta, 0, r, dr, 0, beta, SystemParameters), ThetaSpanEquatorAfter, r_ini_After, dr_ini_After, opts );
% 
%  
% hold on
% [XAfter, YAfter] = pol2cart(ThetaEquatorAfter, rEquatorAfter);
% plot(YAfter, XAfter, 'LineWidth', 3)
% grid on
% axis equal
% % axis([0, 1.5, -0.9, 1])
%         axes = gca;
%         set(gca, 'FontSize', 16)
%         axes.XAxisLocation = 'origin';
%         axes.YAxisLocation = 'origin';
%         

%% Meridional plane: piece-wise integration between the subsolar nose and A1

%%% Meridian Plane: solving explicitly from the subsolar nose to A1 downstream

phiMeridian = pi/2;
dr_ini = 0;
r_ini = rSubSolarNose;
ThetaSpanMeridian = (0 : DeltaThetaDeg : 90-DeltaThetaDeg)*pi/180;

Balance_Meridian_CAN = @(theta, r, drTheta) PressureBalanceCANDisk(theta, phiMeridian, r, drTheta, 0, beta, SystemParameters);

[ThetaMeridian1BeforeA1, rMeridian1BeforeA1] = ode15i(Balance_Meridian_CAN, ThetaSpanMeridian, r_ini, dr_ini );
% 
% 
% hold on
% [X,Y] = pol2cart(ThetaMeridian1BeforeA1, rMeridian1BeforeA1);
% plot(Y, X, 'LineWidth', 3)
% grid on
% axis equal
% % axis([0, 1.3, 0, 1.3])
%         axes = gca;
%         set(gca, 'FontSize', 16)
%         axes.XAxisLocation = 'origin';
%         axes.YAxisLocation = 'origin';


%%% Meridian Plane: solving explicitly from A1 downstream

ThetaSpanMeridianAfterA1 = (90:DeltaThetaDeg:ThetaMaxDeg)*pi/180;
rA1 = rMeridian1BeforeA1(end);

funPreCusp = @(theta, r) drMeridianCAN(theta, r, beta, CAN_DiskParameters);

[thetaMeridian1AfterA1, rMeridian1AfterA1] = ode45( funPreCusp, ThetaSpanMeridianAfterA1, rA1 );
% 
% hold on
% [X,Y] = pol2cart(thetaMeridian1AfterA1, rMeridian1AfterA1);
% plot(Y, X, 'LineWidth', 3)
% grid on
% axis equal
% % axis([0, 1.3, 0, 1.3])
%         axes = gca;
%         set(gca, 'FontSize', 16)
%         axes.XAxisLocation = 'origin';
%         axes.YAxisLocation = 'origin';

%%% Concatenating the two at A1: solving explicitly from the subsolar nose to A1 downstream

rMeridian1 = vertcat(rMeridian1BeforeA1, rMeridian1AfterA1); 
ThetaMeridian1 = vertcat(ThetaMeridian1BeforeA1, thetaMeridian1AfterA1);
% 
% hold on
% [X,Y] = pol2cart(ThetaMeridian1, rMeridian1);
% plot(Y, X, 'LineWidth', 3)
% grid on
% axis equal
% % axis([0, 1.3, 0, 1.3])
%         axes = gca;
%         set(gca, 'FontSize', 16)
%         axes.XAxisLocation = 'origin';
%         axes.YAxisLocation = 'origin';


%% Meridian plane: piece-wise integration at A2

%%% Meridian Plane: finding the intersection r0Y of the surface with the
%%% Y-axis

thetaOY = pi/2;
phiOY = pi/2;
PBMeridian = @(r) 1-2*sqrt(1+beta)*abs( 2*M/r^3 + norm(B_CAN_Meridian(SystemParameters, r, thetaOY, phiOY)) );

rGuess = rSubSolarNose;
r0Y = fzero(PBMeridian, rGuess);


%%% Meridian Plane: solving implicitly from r0Y downstream

phiMeridian = pi/2;
dr_ini_After = 2;
r_ini = r0Y;
ThetaSpanAfter = (90 : DeltaThetaDeg : ThetaMaxDeg)*pi/180;
options = odeset('RelTol',1e-6,'AbsTol',1e-6);

Balance_Meridian_CAN = @(theta, r, drTheta) PressureBalanceCANDisk(theta, phiMeridian, r, drTheta, 0, beta, SystemParameters);

[ThetaAfter, rAfter] = ode15i(Balance_Meridian_CAN, ThetaSpanAfter, r_ini, dr_ini_After, options);

% figure;
% hold on
% [X,Y] = pol2cart(ThetaAfter, rAfter);
% plot(Y, X, 'LineWidth', 3)
% grid on
% axis equal
% %  axis([0, 1.3, 0, 1.3])
%         axes = gca;
%         set(gca, 'FontSize', 16)
%         axes.XAxisLocation = 'origin';
%         axes.YAxisLocation = 'origin';
% 


%%% Meridian Plane: solving implicitly from r0Y upstream
dr_ini_Before = dr_ini_After;
thetaSpanPrevious = (90:-DeltaThetaDeg:DeltaThetaDeg)*pi/180;
[thetaPreviousTemp, rPreviousTemp] = ode15i(Balance_Meridian_CAN, thetaSpanPrevious, r_ini, dr_ini_Before, options);
rPrevious =  vertcat(0, flipud(rPreviousTemp(1:end)));
thetaPrevious = pi/2 - vertcat(thetaPreviousTemp, 0);


% %%% TEST 
% C = 0;
% r_ini_Before = rPrevious(end-C);
% dr_ini_Before = (rPrevious(end-C)-rPrevious(end-C-1))/(DeltaThetaDeg*pi/180);
% ThetaSpanAfter_Before = (90-C*DeltaThetaDeg : DeltaThetaDeg : ThetaMaxDeg)*pi/180;
% [ThetaAfter, rAfter] = ode15i(Balance_Meridian_CAN, ThetaSpanAfter_Before, r_ini_Before, dr_ini_Before, options);

% {
% hold on
% [X,Y] = pol2cart(thetaPrevious, rPrevious);
% plot(Y, X, 'LineWidth', 3)
% grid on
% axis equal
% % axis([0, 1.3, 0, 1.3])
%         axes = gca;
%         set(gca, 'FontSize', 16)
%         axes.XAxisLocation = 'origin';
%         axes.YAxisLocation = 'origin';
% }


%{
% TEST: explicit integrations in 2 directions

%  %%% Meridian Plane:solving explicitely for the POST-cusp section, in 2 directions

% Solving the PB from rOY downstream
% thetaSpanAfter = (90.01:DeltaThetaDeg:ThetaMaxDeg)*(pi/180);
thetaSpanAfter_temp = (90+DeltaThetaDeg:DeltaThetaDeg:ThetaMaxDeg)*(pi/180);
thetaSpanAfter = horzcat(90.01*(pi/180), thetaSpanAfter_temp);

funPostCusp = @(theta, r) drMeridianCAN(theta, r, beta, CAN_DiskParameters);
[thetaAfter, rAfter] = ode45( funPostCusp, thetaSpanAfter, r0Y );
% 

hold on
[X,Y] = pol2cart(thetaAfter, rAfter);
plot(Y, X, 'LineWidth', 3)
grid on
axis equal
%  axis([0, 1.3, 0, 1.3])
        axes = gca;
        set(gca, 'FontSize', 16)
        axes.XAxisLocation = 'origin';
        axes.YAxisLocation = 'origin';

% Solving the PB from rOY upstream
% thetaSpanPrevious = (90-0.001:-DeltaThetaDeg:0)*pi/180;
thetaSpanPrevious_temp = (90-DeltaThetaDeg:-DeltaThetaDeg:0)*pi/180;
thetaSpanPrevious = horzcat((90-0.01)*pi/180, thetaSpanPrevious_temp);

funPreCusp = @(theta, r) drMeridianCAN(theta, r, beta, CAN_DiskParameters);
[thetaPrevious, rPreviousTemp] = ode45( funPreCusp, thetaSpanPrevious, r0Y );

rPrevious = [0; flipud(rPreviousTemp(2:end))];
% Segment upstream of A2


hold on
[X,Y] = pol2cart(thetaPrevious, rPreviousTemp);
plot(Y, X, 'LineWidth', 3)
grid on
axis equal
% axis([0, 1.3, 0, 1.3])
        axes = gca;
        set(gca, 'FontSize', 16)
        axes.XAxisLocation = 'origin';
        axes.YAxisLocation = 'origin';
% 
% hold off
%}

 

%%% Meridian Plane: merging the two curves

ThetaMerged = horzcat( thetaPrevious(1:1:end-1).', ThetaSpanAfter(1:end)).';
% rTwoCurves = vertcat(rPrevious(1:end-1), rAfter);
rTwoCurves = vertcat(rPrevious(1:end), rAfter(2:end));

rMerged =  rTwoCurves;
C=0;
for k = 1:length(ThetaMerged)
    
    if rMerged(k) < rMeridian1(k)
       rMerged(k) = rMeridian1(k);
       C = C+1;
    end
    
end


ThetaCusp = DeltaThetaDeg * ( C-1 ) ;

rMeridian = rMerged;


% figure;
% [X,Y] = pol2cart(ThetaMerged, rMeridian);
% plot(Y, X, 'LineWidth', 3)
% grid on
% axis equal
% % axis([0, 1.3, 0, 1.3])
%         axes = gca;
%         set(gca, 'FontSize', 16)
%         axes.XAxisLocation = 'origin';
%         axes.YAxisLocation = 'origin';



%% Construction of the initial surface

Array = zeros(size(ThetaSpan,2),size(PhiSpan,2));       % Initialising the Array
Array(1,:) = rSubSolarNose * ones(size(PhiSpan));                       % Left-boundary = sub-solar nose
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
%     title('Solution in the equatorial plane')
    hold on
        % Curve
        ThetaAxis = pi/2-ThetaEquator;
        [Yeq, Zeq] = pol2cart(ThetaAxis, rEquator);
        EquatorialSolution = plot(Yeq, Zeq, 'LineWidth', 4);    
        set(EquatorialSolution, 'Color', [00.8500    0.6000    0.4000]); 
        
        % Options
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
        set(Ylabel, 'position', [-0.45 1.35 0]);
        axis equal
        axis([-1 2.5 -1.3 1.5])
        grid on
        box on
        set(gcf,'color','w')
        
        % CAN Disk
        Rp = SystemParameters.Rp;
        r0 = SystemParameters.r0;
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
        DiskColor = [0    0.4500    0.7400];
        bottom = patch(center(1)+[xout,xin], ...
                   center(2)+[yout,yin], ...
                   z1*ones(1,2*length(xout)),'FaceColor', DiskColor, 'FaceAlpha', 0.1, 'LineStyle', '-');
        top = patch(center(1)+[xout,xin], ...
                center(2)+[yout,yin], ...
                z2*ones(1,2*length(xout)),'FaceColor', DiskColor, 'FaceAlpha', 0.1, 'LineStyle', '-');
        [X,Y,Z] = cylinder(1,length(xin));
        outer = surf(rout*X+center(1), ...
                 rout*Y+center(2), ...
                 Z*(z2-z1)+z1,'FaceColor', DiskColor, 'FaceAlpha', 0.5, 'LineStyle', 'none');
        inner = surf(rin*X+center(1), ...
                 rin*Y+center(2), ...
                 Z*(z2-z1)+z1,'FaceColor', DiskColor, 'FaceAlpha', 0.5, 'LineStyle', 'none');
             
        % Annotation: Solar wind and Magnetic moment
        c_SolarWind = [0.6400    0.0800    0.1800];
        c_MagneticMoment = [0.6400    0.0800    0.1800];
        mArrow3([1.30+0.2,0.90+0.2,0],[1.3+0.2,0.7+0.2,0],'color',c_SolarWind,'stemWidth', 0.0085, 'tipWidth', 0.018, 'FaceAlpha', 0.9);
        SolarWind = '$\vec{v}$';
        MagneticMomentVector = '$\vec{\mathcal{M}}$';
        text(1.35+0.2, 0.8+0.2, SolarWind, 'Interpreter','latex', 'FontSize', 20, 'FontName', 'Helvetica', 'color', c_SolarWind);
        rectangle('Position',[-0.06 -0.06 0.12 0.12],'Curvature',[1 1], 'FaceColor', c_MagneticMoment, 'Edgecolor', 'none')
        rectangle('Position',[-0.05 -0.05 0.05*2 0.05*2],'Curvature',[1 1], 'FaceColor', 'white', 'Edgecolor', 'none')
        text(0.1, -0.1, MagneticMomentVector, 'Interpreter','latex', 'FontSize', 16, 'FontWeight', 'bold', 'FontName', 'Helvetica', 'color', c_MagneticMoment);
        line([-0.05, 0.05], [0, 0], 'LineWidth', 1, 'Color', [0.1500    0.1500    0.1500]);
        line([0, 0], [-0.05, 0.05], 'LineWidth', 1, 'Color', [0.1500    0.1500    0.1500]);
        rectangle('Position',[-0.017 -0.017 0.017*2 0.017*2],'Curvature',[1 1], 'FaceColor', c_MagneticMoment, 'Edgecolor', 'none');
        
        % Annotation: Coordinate system
        c_Position = [0.3000    0.4500-0    0.7400];
        line([0,Yeq(100)], [0, Zeq(100)], 'LineWidth', 2, 'LineStyle', ':', 'Color', c_Position);
        circular_arrow(EquatorialSolution, 0.2, [0,0], 65, 30, 1, c_Position, 7, 'vback2', 2)
        PositionTheta = '$\theta$';
        PositionR = '$r$';
        text(0.15-0.07, 0.34-0.07, PositionTheta, 'Interpreter','latex', 'FontSize', 20, 'FontWeight', 'bold', 'FontName', 'Helvetica', 'color', c_Position);
        text(0.53+0.45, 0.35+0.43, PositionR, 'Interpreter','latex', 'FontSize', 22, 'FontWeight', 'bold', 'FontName', 'Helvetica', 'color', c_Position);
	
        % Optional: Add dipolar solution
        load("/Users/flavien/Documents/UCL/Publications/2- GRL - Noon-Midnight Meridian Plane/Figures/DipoleEquatorMeridian.mat");
        [YeqDipole, ZeqDipole] = pol2cart(ThetaAxis, rDipoleEquator);
        EquatorialDipoleSolution_Color = [0.8500    0.6000    0.4000];
        EquatorialDipoleSolution = plot(YeqDipole, ZeqDipole, 'LineWidth', 2, 'LineStyle', '-.', 'Color', EquatorialDipoleSolution_Color);    

    hold off
    
end



%% 
% Piecewise curve solutions in the noon-midnight meridian plane
    
if ismember("CurveMeridianParts", Plots)

    figure;
%     title('Piece-wise solution in the noon-midnight meridian plane')
    hold on
    
        % Curves
        [Y, Z] = pol2cart(pi/2-ThetaMeridian1, rMeridian1);
        MeridianCurve1 = plot(Y, Z, 'LineWidth', 4);
        [Y, Z] = pol2cart(pi/2-ThetaMerged, rTwoCurves);
        MeridianCurve2 = plot(Y, Z, 'LineWidth', 4);
        
        set(MeridianCurve1, 'Color', [0.4000    0.4500    0.7400]); 
        set(MeridianCurve2, 'Color', [0.4700    0.6700    0.5000]); 
        
        % Options
        axes = gca;
        set(axes, 'FontName'   , 'Helvetica',...
                  'XAxisLocation', 'origin',...
                  'YAxisLocation', 'origin',...
                  'FontSize', 16,...
                  'LineWidth'   , 1 ...
              );
        Ylabel = ylabel('Z (r_0)');
        Xlabel = xlabel('Y (r_0)');  
        set([Xlabel, Ylabel], 'FontName', 'AvantGarde', 'FontSize', 18);
        set(Ylabel, 'position', [-0.45 1.35 0]);
        axis equal
        axis([-0.7 2.5 -1.3 1.5])
        grid on
        box on
        set(gcf,'color','w') 
        
        % CAN Disk
        Rp = SystemParameters.Rp;
        r0 = SystemParameters.r0;
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
        DiskColor = [0    0.4500    0.7400];
        bottom = patch(center(1)+[xout,xin], ...
                   center(2)+[yout,yin], ...
                   z1*ones(1,2*length(xout)),'', 'FaceColor', DiskColor, 'FaceAlpha', 0.1, 'LineStyle', '-');
        top = patch(center(1)+[xout,xin], ...
                center(2)+[yout,yin], ...
                z2*ones(1,2*length(xout)),'', 'FaceColor', DiskColor, 'FaceAlpha', 0.1, 'LineStyle', '-');
        [X,Y,Z] = cylinder(1,length(xin));
        outer = surf(rout*X+center(1), ...
                 rout*Y+center(2), ...
                 Z*(z2-z1)+z1, 'FaceColor', DiskColor, 'FaceAlpha', 0.1, 'LineStyle', 'none');
        inner = surf(rin*X+center(1), ...
                 rin*Y+center(2), ...
                 Z*(z2-z1)+z1, 'FaceColor', DiskColor, 'FaceAlpha', 0.1, 'LineStyle', 'none');
        rotate(top, [0 1 0], 90, [0 0 0])
        rotate(bottom, [0 1 0], 90, [0 0 0])
        rotate(outer, [0 1 0], 90, [0 0 0])
        rotate(inner, [0 1 0], 90, [0 0 0])

        % Annotation: Solar wind and Magnetic moment
        c_SolarWind = [0.6400    0.0800    0.1800];
        c_MagneticMoment = [0.6400    0.0800    0.1800];
        mArrow3([1.3+0.2,0.90+0.2,0],[1.3+0.2, 0.7+0.2, 0],'color',c_SolarWind,'stemWidth', 0.0085, 'tipWidth', 0.018, 'FaceAlpha', 0.9);
        SolarWind = '$\vec{v}$';
        MagneticMomentVector = '$\vec{\mathcal{M}}$';
        text(1.35+0.2, 0.8+0.2, SolarWind, 'Interpreter','latex', 'FontSize', 20, 'FontName', 'Helvetica', 'color', c_SolarWind);

        mArrow3([0.18, 0, 0], [-0.18, 0, 0],'color',c_MagneticMoment,'stemWidth', 0.01, 'tipWidth', 0.02);  
        text(0.1, -0.1, MagneticMomentVector, 'Interpreter','latex', 'FontSize', 16, 'FontWeight', 'bold', 'FontName', 'Helvetica', 'color', c_MagneticMoment);

        % Annotation: Coordinate system
        c_Points = [0.7000    0.3300    0.1000];
        A1_PositionInVector = 90/DeltaThetaDeg+1;
        A1_Position = rMeridian1(A1_PositionInVector);
        A1 = plot(A1_Position, 0, 'o', 'color', c_Points, 'MarkerSize', 10, 'LineWidth', 2);
        A2 = plot(r0Y, 0, 'o', 'color', c_Points, 'MarkerSize', 10, 'LineWidth', 2);
        text(A1_Position-0.25, 0.1, '$A_1$', 'Interpreter','latex', 'FontSize', 20, 'FontName', 'Helvetica', 'color', c_Points);
        text(r0Y+0.05, 0.11, '$A_2$', 'Interpreter','latex', 'FontSize', 20, 'FontName', 'Helvetica', 'color', c_Points);
        
        
    hold off
    
end



%%
% Curve solution in the noon-midnight meridian plane

if ismember("CurveMeridian", Plots)
    
    figure;
    title('Solution in the noon-midnight meridian plane')
    hold on
    
        % Curve
        [Y, Z] = pol2cart(pi/2-ThetaMerged, rMerged);
        MeridianSolution = plot(Y, Z, 'LineWidth', 4);
        
        PointXCusp = [0, (rSubSolarNose-0.05)*sin( ThetaCusp*pi/180 )];
        PointYCusp = [0  (rSubSolarNose-0.05)*cos( ThetaCusp*pi/180 )]; 
        CuspPosition = plot(PointXCusp, PointYCusp, ':', 'LineWidth', 2);
        
        MeridianSolutionColor = [0.4000    0.6000    0.3000];
        set(MeridianSolution, 'Color', MeridianSolutionColor); 
        set(CuspPosition, 'Color', [0.4900    0.1800    0.5000]); 
        
        % Options
        axes = gca;
        set(axes, 'FontName'   , 'Helvetica',...
                  'XAxisLocation', 'origin',...
                  'YAxisLocation', 'origin',...
                  'FontSize', 16,...
                  'LineWidth'   , 1 ...
              );
        Ylabel = ylabel('Z (r_0)');
        Xlabel = xlabel('Y (r_0)');  
        set([Xlabel, Ylabel], 'FontName', 'AvantGarde', 'FontSize', 18);
        set(Ylabel, 'position', [-0.45 1.35 0]);
        axis equal
         axis([-0.7 2.5 -1.3 1.5])
        grid on
        box on
        set(gcf,'color','w') 
        

        % CAN Disk
        Rp = SystemParameters.Rp;
        r0 = SystemParameters.r0;
        
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
        DiskColor = [0    0.4500    0.7400];
        bottom = patch(center(1)+[xout,xin], ...
                   center(2)+[yout,yin], ...
                   z1*ones(1,2*length(xout)),'', 'FaceColor', DiskColor, 'FaceAlpha', 0.1, 'LineStyle', '-');
        top = patch(center(1)+[xout,xin], ...
                center(2)+[yout,yin], ...
                z2*ones(1,2*length(xout)),'', 'FaceColor', DiskColor, 'FaceAlpha', 0.1, 'LineStyle', '-');
        [X,Y,Z] = cylinder(1,length(xin));
        outer = surf(rout*X+center(1), ...
                 rout*Y+center(2), ...
                 Z*(z2-z1)+z1, 'FaceColor', DiskColor, 'FaceAlpha', 0.1, 'LineStyle', 'none');
        inner = surf(rin*X+center(1), ...
                 rin*Y+center(2), ...
                 Z*(z2-z1)+z1, 'FaceColor', DiskColor, 'FaceAlpha', 0.1, 'LineStyle', 'none');
        rotate(top, [0 1 0], 90, [0 0 0])
        rotate(bottom, [0 1 0], 90, [0 0 0])
        rotate(outer, [0 1 0], 90, [0 0 0])
        rotate(inner, [0 1 0], 90, [0 0 0])

        % Annotation: Solar wind and Magnetic moment
        c_SolarWind = [0.6400    0.0800    0.1800];
        c_MagneticMoment = [0.6400    0.0800    0.1800];
        mArrow3([1.3+0.2,0.90+0.2,0],[1.3+0.2, 0.7+0.2, 0],'color',c_SolarWind,'stemWidth', 0.0085, 'tipWidth', 0.018, 'FaceAlpha', 0.9);
        SolarWind = '$\vec{v}$';
        MagneticMomentVector = '$\vec{\mathcal{M}}$';
        text(1.35+0.2, 0.8+0.2, SolarWind, 'Interpreter','latex', 'FontSize', 20, 'FontName', 'Helvetica', 'color', c_SolarWind);

        mArrow3([0.18, 0, 0], [-0.18, 0, 0],'color',c_MagneticMoment,'stemWidth', 0.01, 'tipWidth', 0.02);  
        text(0.1, -0.1, MagneticMomentVector, 'Interpreter','latex', 'FontSize', 16, 'FontWeight', 'bold', 'FontName', 'Helvetica', 'color', c_MagneticMoment);
        
        % Annotation: Coordinate system
        c_Points = [0.7000    0.3300    0.1000];
        c_Position = [0.4900    0.1800    0.5000];
        circular_arrow(MeridianSolution, 0.3, [0,0], 55, 50, 1, c_Position, 7, 'vback2', 2);
        PositionTheta = '$\theta_C$';
        text(0.2, 0.32, PositionTheta, 'Interpreter','latex', 'FontSize', 20, 'FontWeight', 'bold', 'FontName', 'Helvetica', 'color', c_Position);
        A2 = plot(r0Y, 0, 'o', 'color', c_Points, 'MarkerSize', 10, 'LineWidth', 2);
        text(r0Y+0.05, 0.11, '$A_2$', 'Interpreter','latex', 'FontSize', 20, 'FontName', 'Helvetica', 'color', c_Points);
        text(PointXCusp(2)+0.05, PointYCusp(2)+0.05, '$C$', 'Interpreter','latex', 'FontSize', 20, 'FontName', 'Helvetica', 'color', c_Position);
        
        % Optional: Add dipolar solution
        load("/Users/flavien/Documents/UCL/Publications/2- GRL - Noon-Midnight Meridian Plane/Figures/DipoleEquatorMeridian.mat");
        ThetaAxis = pi/2-ThetaEquator;
        [YMeridianDipole, ZeqMeridianDipole] = pol2cart(ThetaAxis, rDipoleMeridian);
        MeridianDipoleSolution_Color = horzcat(MeridianSolutionColor, 0.3);
        MeridianDipoleSolution = plot(YMeridianDipole, ZeqMeridianDipole, 'LineWidth', 2, 'LineStyle', '-.', 'Color', MeridianDipoleSolution_Color);

        
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
             Z*(z2-z1)+z1, 'FaceAlpha', 0.5, 'LineStyle', 'none');
    inner = surf(rin*X+center(1), ...
             rin*Y+center(2), ...
             Z*(z2-z1)+z1, 'FaceAlpha', 0.5, 'LineStyle', 'none');

    
    ConstructionSurface(Array,ThetaSpan,PhiSpan)
    set(gcf,'color','w');
    axis([-1.5, 1.5, -2.3, 2.3, -2.3, 2.3])
    box on 
    view([45+90+5 30])
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
         axis([-2, 2, -2, 2, -2, 2]);
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
               z1*ones(1,2*length(xout)),'', 'FaceAlpha', 0.5);
    top = patch(center(1)+[xout,xin], ...
            center(2)+[yout,yin], ...
            z2*ones(1,2*length(xout)),'', 'FaceAlpha', 0.5);
    [X,Y,Z] = cylinder(1,length(xin));
    outer = surf(rout*X+center(1), ...
             rout*Y+center(2), ...
             Z*(z2-z1)+z1, 'FaceAlpha', 0.5);
    inner = surf(rin*X+center(1), ...
             rin*Y+center(2), ...
             Z*(z2-z1)+z1, 'FaceAlpha', 0.5);
    rotate(inner, [1, 0, 0], 90 )
    rotate(outer, [1, 0, 0], 90 )
    rotate(bottom, [1, 0, 0], 90 )
    rotate(top, [1, 0, 0], 90 )
    
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
    FValuesGridInsideVector = F_values(rGuess, rEquator, rMeridian, rSubSolarNose, ThetaMaxDeg, PhiMaxDeg, DeltaThetaDeg, DeltaPhiDeg, NbPointsGridInside, beta, SystemParameters);
    FValuesGridInside = reshape(FValuesGridInsideVector,[NbPointsPhiInside, NbPointsThetaInside]);
    ErrorAbsolute = log10(abs(FValuesGridInside(1:1:end,1:1:end)));
    pcolor(ErrorAbsolute(1:end, 1:end)),shading flat; colorbar; 
    caxis([-4 , 0])
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
    FValuesGridInsideVector = F_values(rGuess, rEquator, rMeridian, rSubSolarNose, ThetaMaxDeg, PhiMaxDeg, DeltaThetaDeg, DeltaPhiDeg, NbPointsGridInside, beta, SystemParameters);
    FValuesGridInside = reshape(FValuesGridInsideVector,[NbPointsPhiInside, NbPointsThetaInside]);
    ErrorAbsolute = log10(abs(FValuesGridInside(1:1:end,1:1:end)));
    pcolor(ErrorAbsolute(1:end, 1:end));shading flat;
    MaxErrorScale = log10(2);
    MinErrorScale = -4;
        caxis([MinErrorScale MaxErrorScale]);
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
%     zz = 3*Rp/r0*ones(size(xx));
    surf(yy,zz,xx, 'FaceAlpha', 0.3)
    
    ImgRGB = imread('ColourMap.png');
    [ImRows, ImCols, ImPlanes] = size(ImgRGB);
    ConstructionSurfaceWarp(Array(1:1:end, 1:1:end),ThetaSpan(1:1:end), PhiSpan(1:1:end), ImgRGB)
    set(gcf,'color','w');
%     colorbar('southoutside');
    colorbar('southoutside', 'YTickLabel', {'0.01', '', '0.1', '', '1', '', '10', '', '100'});
    caxis([MinErrorScale MaxErrorScale]);
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
    surf(yy,zz,xx, 'FaceAlpha', 0.3)
    
    ImgRGB = imread('ColourMap.png');
    [ImRows, ImCols, ImPlanes] = size(ImgRGB);
    ConstructionSurfaceWarp(Array(1:1:end, 1:1:end),ThetaSpan(1:1:end), PhiSpan(1:1:end), ImgRGB)
    set(gcf,'color','w');
%     colorbar('southoutside');
    colorbar('southoutside', 'YTickLabel', {'0.01', '', '0.1', '', '1', '', '10', '', '100'});
        caxis([-4 0]);
    view([0 90])
    
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
    title('Relative PB assessment projected onto the initial surface')  
    ImgRGB = imread('ColourMap.png');
    [ImRows, ImCols, ImPlanes] = size(ImgRGB);
    ConstructionSurfaceWarp(Array(1:1:end, 1:1:end),ThetaSpan(1:1:end), PhiSpan(1:1:end), ImgRGB)
    set(gcf,'color','w');
%     colorbar('southoutside');
    colorbar('southoutside', 'YTickLabel', {'0.01', '', '0.1', '', '1', '', '10', '', '100'});
    caxis([-4 0]);
    view([0 90])
    axis([-1.7, 1.7, -1.7, 1.7, -1.7, 1.7])
    box on
    grid on

end


end