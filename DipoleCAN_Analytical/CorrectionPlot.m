function CorrectionPlot = CorrectionPlot(Plots, DeltaThetaDeg, ThetaCusp, DeltaPhiDeg, ThetaMaxDeg, PhiMaxDeg, rEquator, rMeridian, SurfaceCorrected, rCroppedCorrected, SurfaceCroppedExtrapolated)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs
%    - Plots: list of strings, selection of plots to display
%           - "GridCroppedAnalysis": colourmap assessing the pressure balance on
%               the planar upper sub-grid
%           - "SurfaceCroppedCorrected: corrected 3D cropped surface on the upper sub-grid
%           - "SurfaceTotCorrected": corrected 3D total surface
%           - "GridWrapped": colourmap assessing the pressure balance
%                   projected onto the corrected 3D surface
%    - DeltaThetaDeg, DeltaPhiDeg: scalars, angular increments, in degrees
%    - ThetaCusp: scalar, theta value of the cusp, in degrees
%    - ThetaMaxDeg, PhiMaxDeg: scalars, maximum values for theta and phi,
%            in degrees
%    - rEquator: vector, solution in the equatorial plane
%           theta = 0:DeltaThetaDeg:ThetaMaxDeg
%    - rMeridian: vector, solution in the noon-midnight meridian plane
%           theta = 0:DeltaThetaDeg:ThetaMaxDeg
%    - SurfaceCorrected: phi*theta array, corrected surface
%           on the entire grid
%     - rCroppedCorrected: vector, corrected upper sub-grid
%     - SurfaceCroppedExtrapolated: phi*theta array, corrected surface on
%     the sub-upper grid
%
% Output
%     - Graphes, as indicated by Plots 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%
% Absolute pressure balance assessment on the cropped 2D grid

if ismember("GridCroppedAnalysis", Plots)
    
    ThetaMaxDegSubSolar =  (ThetaCusp -  mod(ThetaCusp, DeltaThetaDeg)) - 0*DeltaThetaDeg;
    PhiMaxDegSubSolar = PhiMaxDeg;
    NbPointsThetaSubSolar = (ThetaMaxDegSubSolar/DeltaThetaDeg+1);                         % Number of points in theta-direction for phi=cst (theta = 0:Max_deg_equ)
    NbPointsPhiSubSolar = (PhiMaxDeg/DeltaPhiDeg+1);                                       % Number of points in phi-direction for theta=cst (phi = 0:Max_deg_phi)
    NbPointsGridSubSolarTotal = NbPointsThetaSubSolar*NbPointsPhiSubSolar;
    NbPointsThetaSubSolarInside = NbPointsThetaSubSolar-2;                               % Number of points on the equator without first (theta=0) and last (theta=theta_max) points
    NbPointsPhiSubSolarInside = NbPointsPhiSubSolar-2;
    NbPointsGridSubSolarInside = NbPointsThetaSubSolarInside*NbPointsPhiSubSolarInside;  % Number of points on the interior grid for theta, phi = delta : Max_deg-delta
       
    figure;
    FValuesGridInsideVectorSubSolar = F_valuesRelative(rCroppedCorrected, rEquator, rMeridian, 1, ThetaMaxDegSubSolar, PhiMaxDegSubSolar, DeltaThetaDeg, DeltaPhiDeg, NbPointsGridSubSolarInside);
    FValuesGridInsideSubSolar=reshape(FValuesGridInsideVectorSubSolar,[NbPointsPhiSubSolarInside, NbPointsThetaSubSolarInside]);
    ErrorAbsoluteSubSolar = (abs(FValuesGridInsideSubSolar(1:1:end,1:1:end)));
    pcolor((ErrorAbsoluteSubSolar(1:end, 1:end))),shading flat; colorbar; 
    title('Absolute PB values on the corrected cropped surface')  
%     caxis([0 0.13]);
    axes = gca;
    xlabel('Theta','FontSize',16,'FontWeight','bold');
    xticks(0:10/DeltaThetaDeg:NbPointsThetaSubSolarInside);
    xticklabels({(0:10/DeltaThetaDeg:NbPointsThetaSubSolarInside)*DeltaThetaDeg});
    ylabel('Phi','FontSize',16,'FontWeight','bold');
    yticks(0:10/DeltaPhiDeg:NbPointsPhiSubSolarInside);
    yticklabels({(0:10/DeltaPhiDeg:NbPointsPhiSubSolarInside)*DeltaPhiDeg});
    set(gca, 'FontSize', 16);
    axes.XAxisLocation = 'origin';
    axes.YAxisLocation = 'origin';
    set(gcf,'color','w');
    box on
    
end



%%
% Corrected cropped 2D surface

if ismember("SurfaceCroppedCorrected", Plots)
    
    ThetaMaxDegSubSolar =  (ThetaCusp -  mod(ThetaCusp, DeltaThetaDeg)) - 0*DeltaThetaDeg;
    PhiMaxDegSubSolar = PhiMaxDeg;
    Offset = (ThetaCusp - ThetaMaxDegSubSolar + DeltaThetaDeg)/DeltaThetaDeg;    
    NbPointsThetaSubSolar = (ThetaMaxDegSubSolar/DeltaThetaDeg+1);                         % Number of points in theta-direction for phi=cst (theta = 0:Max_deg_equ)
    NbPointsPhiSubSolar = (PhiMaxDeg/DeltaPhiDeg+1);                                       % Number of points in phi-direction for theta=cst (phi = 0:Max_deg_phi)
    NbPointsGridSubSolarTotal = NbPointsThetaSubSolar*NbPointsPhiSubSolar;
    NbPointsThetaSubSolarInside = NbPointsThetaSubSolar-2;                               % Number of points on the equator without first (theta=0) and last (theta=theta_max) points
    NbPointsPhiSubSolarInside = NbPointsPhiSubSolar-2;
    NbPointsGridSubSolarInside = NbPointsThetaSubSolarInside*NbPointsPhiSubSolarInside;  % Number of points on the interior grid for theta, phi = delta : Max_deg-delta

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
    surf(yy,zz,xx, 'FaceAlpha', 0.3)

    % Surface
    title('Corrected cropped surface')  
    ThetaSpanSubSolarExtrap = (0:DeltaThetaDeg:ThetaMaxDegSubSolar-DeltaThetaDeg+Offset*DeltaThetaDeg)*pi/180;
    PhiSpanSubSolarExtrap = (0:DeltaPhiDeg:PhiMaxDegSubSolar)*pi/180;
    ConstructionSurface(SurfaceCroppedExtrapolated.',ThetaSpanSubSolarExtrap, PhiSpanSubSolarExtrap)
    axis equal
    set(gcf,'color','w');
    axis([-1.2, 1.2, -1.2, 1.2, -0.2, 1.2])

end



%%
% Corrected total 2D surface

if ismember("SurfaceTotCorrected", Plots)
    
    figure;
    title('Corrected surface')
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
    
%     surf(yy,zz,xx, 'FaceAlpha', 0.3)
    surf(xx, yy,zz, 'FaceAlpha', 0.8)
    PhiSpanConcatenated = (0:DeltaPhiDeg:PhiMaxDeg)*pi/180;
    ThetaSpanConcatenated = (0:DeltaThetaDeg:ThetaMaxDeg)*pi/180;
    ConstructionSurface(SurfaceCorrected(1:1:end, 1:1:end).',ThetaSpanConcatenated(1:1:end), PhiSpanConcatenated(1:1:end))
   box on
   axis([-1.3 1.3 -1.7 1.7 -1.7 1.7])
%    view(45+90, 45)
   grid on

end



%%
% Iso-theta contours of the corrected total surface 

if ismember("SurfaceTotCorrectedContours", Plots)
    
    figure;
    title('Corrected surface contours')
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
    
    PhiSpanConcatenated = (0:DeltaPhiDeg:PhiMaxDeg)*pi/180;
    ThetaSpanConcatenated = (0:DeltaThetaDeg:ThetaMaxDeg)*pi/180;
    ConstructionSurfaceContours(SurfaceCorrected(1:1:end, 1:1:end).',ThetaSpanConcatenated(1:1:end), PhiSpanConcatenated(1:1:end))
   box on
   grid on

end



%%
% Corrected total surface + Field lines

if ismember("SurfaceTotCorrectedFieldLines", Plots)

    figure;
    title('Corrected surface and dipole field lines')  

    L = (0.1:0.5:50).';
    theta = (0:0.1:180)*pi/180;
    phi = 90*pi/180;

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
    PhiSpanConcatenated = (0:DeltaPhiDeg:PhiMaxDeg)*pi/180;
    ThetaSpanConcatenated = (0:DeltaThetaDeg:ThetaMaxDeg)*pi/180;
    ConstructionSurface(SurfaceCorrected(1:1:end, 1:1:end).',ThetaSpanConcatenated(1:1:end), PhiSpanConcatenated(1:1:end))
    
end
    
    

%%
% Absolute pressure balance assessment projected onto the corrected surface

if ismember("GridWrappedAbsolute", Plots)
    
    NbPointsTheta = ThetaMaxDeg/DeltaThetaDeg+1;   % Number of points in theta-direction for phi=cst (theta = 0:Max_deg_equ)
    NbPointsPhi = PhiMaxDeg/DeltaPhiDeg+1;         % Number of points in phi-direction for theta=cst (phi = 0:Max_deg_phi)
    NbPointsGrid = NbPointsTheta*NbPointsPhi;
    NbPointsThetaInside = NbPointsTheta-2;                           % Number of points on the equator without first (theta=0) and last (theta=theta_max) points
    NbPointsPhiInside = NbPointsPhi-2;
    NbPointsGridInside = NbPointsThetaInside*NbPointsPhiInside;      % Number of points on the interior grid for theta, phi = delta : Max_deg-delta

    PhiSpanConcatenated = (0:DeltaPhiDeg:PhiMaxDeg)*pi/180;
    ThetaSpanConcatenated = (0:DeltaThetaDeg:ThetaMaxDeg)*pi/180;
    
    figure();
    CorrectedGridInside = SurfaceCorrected(2:end-1, 2:end-1);
    ConcatenatedSurfaceInterpVec = CorrectedGridInside(:);
    FValuesConcatenatedGridInsideVector = F_values(ConcatenatedSurfaceInterpVec, rEquator, rMeridian, 1, ThetaMaxDeg, PhiMaxDeg, DeltaThetaDeg, DeltaPhiDeg, NbPointsGridInside, CAN_DiskParameters);
    FValuesConcatenatedGridInside=reshape(FValuesConcatenatedGridInsideVector,[NbPointsPhiInside, NbPointsThetaInside]);
    ErrorAbsolute = (abs(FValuesConcatenatedGridInside(1:1:end,1:1:end)));
    pcolor(ErrorAbsolute(1:end, 1:end)),shading flat
    MaxErrorScale = 0.03;
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
    close()
        
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
    surf(yy,zz,xx, 'FaceAlpha', 0.3)
    
    title('Absolute PB assessment projected onto the corrected surface')  
    
    ImgRGB = imread('ColourMap.png');
    [ImRows, ImCols, ImPlanes] = size(ImgRGB);
    ConstructionSurfaceWarp(SurfaceCorrected(1:1:end, 1:1:end).',ThetaSpanConcatenated(1:1:end), PhiSpanConcatenated(1:1:end), ImgRGB)
    set(gcf,'color','w');
%     colorbar('southoutside');
    colorbar('southoutside', 'YTickLabel', {'0.01', '', '0.1', '', '1', '', '10', '', '100'});
    caxis([0 MaxErrorScale]);
    view([0 90])   
    
end



%%
% Absolute pressure balance assessment projected onto the corrected surface
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
    CorrectedGridInside = SurfaceCorrected(2:end-1, 2:end-1);
    ConcatenatedSurfaceInterpVec = CorrectedGridInside(:);
    FValuesConcatenatedGridInsideVector = F_values(ConcatenatedSurfaceInterpVec, rEquator, rMeridian, 1, ThetaMaxDeg, PhiMaxDeg, DeltaThetaDeg, DeltaPhiDeg, NbPointsGridInside);
    FValuesConcatenatedGridInside=reshape(FValuesConcatenatedGridInsideVector,[NbPointsPhiInside, NbPointsThetaInside]);
    ErrorAbsolute = (abs(FValuesConcatenatedGridInside(1:1:end,1:1:end)));
    pcolor(ErrorAbsolute(1:end, 1:end)),shading flat
    MaxErrorScale = 0.03;
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
    close()
        
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
    title('Absolute PB and magnetic field lines') 
    ImgRGB = imread('ColourMap.png');
    [ImRows, ImCols, ImPlanes] = size(ImgRGB);
    ConstructionSurfaceWarp(SurfaceCorrected(1:1:end, 1:1:end).',ThetaSpanConcatenated(1:1:end), PhiSpanConcatenated(1:1:end), ImgRGB)
    
    set(gcf,'color','w');
%     colorbar('southoutside');
    colorbar('southoutside', 'YTickLabel', {'0.01', '', '0.1', '', '1', '', '10', '', '100'});
    caxis([0 MaxErrorScale]);
    view([0 90]) 
    axis([-1.7, 1.7, -1.7, 1.7, -1.7, 1.7])

end



%%
% Relative pressure balance assessment projected onto the corrected surface

if ismember("GridWrappedRelative", Plots)
    
    NbPointsTheta = ThetaMaxDeg/DeltaThetaDeg+1;   % Number of points in theta-direction for phi=cst (theta = 0:Max_deg_equ)
    NbPointsPhi = PhiMaxDeg/DeltaPhiDeg+1;         % Number of points in phi-direction for theta=cst (phi = 0:Max_deg_phi)
    NbPointsGrid = NbPointsTheta*NbPointsPhi;
    NbPointsThetaInside = NbPointsTheta-2;                           % Number of points on the equator without first (theta=0) and last (theta=theta_max) points
    NbPointsPhiInside = NbPointsPhi-2;
    NbPointsGridInside = NbPointsThetaInside*NbPointsPhiInside;      % Number of points on the interior grid for theta, phi = delta : Max_deg-delta

    PhiSpanConcatenated = (0:DeltaPhiDeg:PhiMaxDeg)*pi/180;
    ThetaSpanConcatenated = (0:DeltaThetaDeg:ThetaMaxDeg)*pi/180;
    
    figure();
    CorrectedGridInside = SurfaceCorrected(2:end-1, 2:end-1);
    ConcatenatedSurfaceInterpVec = CorrectedGridInside(:);
    FValuesConcatenatedGridInsideVector = F_valuesRelative(ConcatenatedSurfaceInterpVec, rEquator, rMeridian, 1, ThetaMaxDeg, PhiMaxDeg, DeltaThetaDeg, DeltaPhiDeg, NbPointsGridInside);
    FValuesConcatenatedGridInside=reshape(FValuesConcatenatedGridInsideVector,[NbPointsPhiInside, NbPointsThetaInside]);
    ErrorAbsolute = log10(abs(FValuesConcatenatedGridInside(1:1:end,1:1:end)));
    pcolor(ErrorAbsolute(1:end, 1:end)); shading flat;
    MaxErrorScale = 0.05*0;
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
    close()
        
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
    surf(yy,zz,xx, 'FaceAlpha', 0.3)
    
    title('Relative PB assessment projected onto the corrected surface')  
    
    ImgRGB = imread('ColourMap.png');
    [ImRows, ImCols, ImPlanes] = size(ImgRGB);
    ConstructionSurfaceWarp(SurfaceCorrected(1:1:end, 1:1:end).',ThetaSpanConcatenated(1:1:end), PhiSpanConcatenated(1:1:end), ImgRGB)
    set(gcf,'color','w');
    colorbar('southoutside', 'YTickLabel', {'0.01', '', '0.1', '', '1', '', '10', '', '100'});
%     colorbar('southoutside');
    caxis([-4 0]);
    view([0 90])   
    
    % TEST
%     ConstructionSurfaceContours(SurfaceCorrected(1:1:end, 1:1:end).',ThetaSpanConcatenated(1:1:end), PhiSpanConcatenated(1:1:end))

end



%%
% Absolute pressure balance assessment projected onto the corrected surface
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
    CorrectedGridInside = SurfaceCorrected(2:end-1, 2:end-1);
    ConcatenatedSurfaceInterpVec = CorrectedGridInside(:);
    FValuesConcatenatedGridInsideVector = F_valuesRelative(ConcatenatedSurfaceInterpVec, rEquator, rMeridian, 1, ThetaMaxDeg, PhiMaxDeg, DeltaThetaDeg, DeltaPhiDeg, NbPointsGridInside);
    FValuesConcatenatedGridInside=reshape(FValuesConcatenatedGridInsideVector,[NbPointsPhiInside, NbPointsThetaInside]);
    ErrorAbsolute = log10(abs(FValuesConcatenatedGridInside(1:1:end,1:1:end)));
    pcolor(ErrorAbsolute(1:end, 1:end)),shading flat
    MaxErrorScale = 0.05;
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
    close()
        
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
    title('Relative PB and magnetic field lines') 
    ImgRGB = imread('ColourMap.png');
    [ImRows, ImCols, ImPlanes] = size(ImgRGB);
    ConstructionSurfaceWarp(SurfaceCorrected(1:1:end, 1:1:end).',ThetaSpanConcatenated(1:1:end), PhiSpanConcatenated(1:1:end), ImgRGB)
    
    set(gcf,'color','w');
%     colorbar('southoutside');
    colorbar('southoutside', 'YTickLabel', {'0.01', '', '0.1', '', '1', '', '10', '', '100'});
    caxis([-4 0]);
    view([0 90]) 
    axis([-1.7, 1.7, -1.7, 1.7, -1.7, 1.7])

end



end