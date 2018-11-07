function CorrectionPlot = CorrectionPlot(Plots, DeltaThetaDeg, DeltaPhiDeg, ThetaMaxDeg, PhiMaxDeg, rSubSolarNose, SurfaceCorrected, SystemParameters)


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
% Corrected total 2D surface

if ismember("SurfaceCorrected", Plots)
    
    figure;
    title('Corrected surface')
    hold on
    
        % Ring
        Rp = SystemParameters.Rp;        % Planet radius, m
        r0 = SystemParameters.r0;          % Distance scale
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
        a = (SystemParameters.a)*Rp/r0;      % Inner radius of the current distribution (in r0)
        b = (SystemParameters.b)*Rp/r0;   % Outer radius of the current distribution (in r0)
        D = (SystemParameters.D)*Rp/r0;       % Disk Semi-thickness (in r0)
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


        PhiSpanConcatenated = (0:DeltaPhiDeg:PhiMaxDeg)*pi/180;
        ThetaSpanConcatenated = (0:DeltaThetaDeg:ThetaMaxDeg)*pi/180;
        Delta = 2;
        ConstructionSurfaceContours(SurfaceCorrected(1:1:end, 1:1:end).',ThetaSpanConcatenated(1:1:end), PhiSpanConcatenated(1:1:end), Delta)
       box on
       axis([-2.3 2.3 -2.3 2.3 -2.3 2.3])
       view(45+90, 30)
       grid on
       
	hold off

end


%%
% Corrected total surface + Field lines

if ismember("SurfaceCorrectedFieldLines", Plots)

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
             plot3(Z, X, Y, ':', 'LineWidth', 2)
          end
        end
    
    
        % Ring
        Rp = SystemParameters.Rp;        % Planet radius, m
        r0 = SystemParameters.r0;          % Distance scale
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
        a = (SystemParameters.a)*Rp/r0;      % Inner radius of the current distribution (in r0)
        b = (SystemParameters.b)*Rp/r0;   % Outer radius of the current distribution (in r0)
        D = (SystemParameters.D)*Rp/r0;       % Disk Semi-thickness (in r0)
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


        PhiSpanConcatenated = (0:DeltaPhiDeg:PhiMaxDeg)*pi/180;
        ThetaSpanConcatenated = (0:DeltaThetaDeg:ThetaMaxDeg)*pi/180;
        Delta = 2;
        ConstructionSurfaceContours(SurfaceCorrected(1:1:end, 1:1:end).',ThetaSpanConcatenated(1:1:end), PhiSpanConcatenated(1:1:end), Delta)
       box on
       axis([-2.3 2.3 -2.3 2.3 -2.3 2.3])
       view(45+90, 30)
       grid on
       
	hold off
    
end
    
    

%%
% Absolute pressure balance assessment projected onto the corrected surface

if ismember("GridWrapped", Plots)
    
    Jump = 2;
    
    [NbPointsGrid, NbPointsGridInside,      ...
    NbPointsThetaInside, NbPointsPhiInside, ...
    ThetaSpanVector, PhiSpanVector]         = GridDetails(ThetaMaxDeg, PhiMaxDeg, DeltaThetaDeg, DeltaPhiDeg);
    
    figure();
    SurfaceCorrectedInside = SurfaceCorrected(2:end-1, 2:end-1);
    rTotalSurface = SurfaceCorrectedInside(:);
    rBottom = SurfaceCorrected(1,:);
    rTop = SurfaceCorrected(end,:);
    figure;
    PB_Assessment = F_values_Numerical(rTotalSurface, rBottom, rTop, rSubSolarNose, ThetaMaxDeg, PhiMaxDeg, DeltaThetaDeg, DeltaPhiDeg, NbPointsGridInside, SystemParameters);
    ErrorVector = log10(abs(PB_Assessment));
    ErrorGrid = reshape(ErrorVector, NbPointsPhiInside, NbPointsThetaInside);
    pcolor(ErrorGrid(:, Jump:Jump:end));
    shading flat
        caxis([-4 0])
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
        
    hold on

    % Ring
        Rp = SystemParameters.Rp;        % Planet radius, m
        r0 = SystemParameters.r0;          % Distance scale
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
        a = (SystemParameters.a)*Rp/r0;      % Inner radius of the current distribution (in r0)
        b = (SystemParameters.b)*Rp/r0;   % Outer radius of the current distribution (in r0)
        D = (SystemParameters.D)*Rp/r0;       % Disk Semi-thickness (in r0)
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
        Opacity = 0.2;
        bottom = patch(center(1)+[xout,xin], ...
                   center(2)+[yout,yin], ...
                   z1*ones(1,2*length(xout)),'', 'FaceAlpha', Opacity, 'LineStyle', 'none');
        top = patch(center(1)+[xout,xin], ...
                center(2)+[yout,yin], ...
                z2*ones(1,2*length(xout)),'', 'FaceAlpha', Opacity, 'LineStyle', 'none');
        [X,Y,Z] = cylinder(1,length(xin));
        outer = surf(rout*X+center(1), ...
                 rout*Y+center(2), ...
                 Z*(z2-z1)+z1, 'FaceColor', 'blue', 'FaceAlpha', 0.1, 'LineStyle', 'none');
        inner = surf(rin*X+center(1), ...
                 rin*Y+center(2), ...
                 Z*(z2-z1)+z1, 'FaceColor', 'blue', 'FaceAlpha', 0.1, 'LineStyle', 'none');

        PhiSpanConcatenated = (0:DeltaPhiDeg:PhiMaxDeg)*pi/180;
        ThetaSpanConcatenated = (0:DeltaThetaDeg:ThetaMaxDeg)*pi/180;

            ImgRGB = imread('ColourMap.png');
            [ImRows, ImCols, ImPlanes] = size(ImgRGB);
            ConstructionSurfaceWarp(SurfaceCorrected(:,Jump:Jump:end).', ThetaSpanConcatenated(Jump:Jump:end), PhiSpanConcatenated(1:1:end), ImgRGB)
            set(gcf,'color','w');
            view([45+3*90+5 30])
            axis([-2.3, 2.3, -2.3, 2.3, -2.3, 2.3])
            box on
            grid on

        caxis([-4 0])
        cb = colorbar('eastoutside', 'YTickLabel', {'0.01', '', '0.1', '', '1', '', '10', '', '100'}, 'LineWidth', 1, 'FontName' , 'Helvetica', 'FontSize', 15);
        
        % Optional colorbar positioning for publications
%         cbPosition = cb.Position;
%         ZStretch = 0.8;
%         set(cb, 'Position', [cbPosition(1)*0.99 cbPosition(2)+(cbPosition(4)*(1-ZStretch))/2 cbPosition(3) cbPosition(4)*ZStretch])
            
	hold off
    
end



%%
% Absolute pressure balance assessment projected onto the corrected surface
% + Field lines

if ismember("GridWrappedFieldLines", Plots)
    
    Jump = 2;
    
    [NbPointsGrid, NbPointsGridInside,      ...
    NbPointsThetaInside, NbPointsPhiInside, ...
    ThetaSpanVector, PhiSpanVector]         = GridDetails(ThetaMaxDeg, PhiMaxDeg, DeltaThetaDeg, DeltaPhiDeg);
    
    figure();
    SurfaceCorrectedInside = SurfaceCorrected(2:end-1, 2:end-1);
    rTotalSurface = SurfaceCorrectedInside(:);
    rBottom = SurfaceCorrected(1,:);
    rTop = SurfaceCorrected(end,:);
    figure;
    PB_Assessment = F_values_Numerical(rTotalSurface, rBottom, rTop, rSubSolarNose, ThetaMaxDeg, PhiMaxDeg, DeltaThetaDeg, DeltaPhiDeg, NbPointsGridInside, SystemParameters);
    ErrorVector = log10(abs(PB_Assessment));
    ErrorGrid = reshape(ErrorVector, NbPointsPhiInside, NbPointsThetaInside);
    pcolor(ErrorGrid(:, Jump:Jump:end));
    shading flat
        caxis([-4 0])
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
        
    hold on
    
    % L-Shells
     L = (0.1:0.5:20).';
    theta = (0:0.1:180)*pi/180;
    phi = 90*pi/180;
    
    for m = 1:length(phi)
      for k = 1:length(L)
         r = L(k).*sin(theta).^2;
         [X, Y, Z] = sph2cart(phi(m), theta, r);
         axis([-1.7, 1.7, -1.7, 1.7, -1.7, 1.7]);
         plot3(Z, X, Y, ':', 'LineWidth', 2)
      end
    end

    % Ring
        Rp = SystemParameters.Rp;        % Planet radius, m
        r0 = SystemParameters.r0;          % Distance scale
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
        a = (SystemParameters.a)*Rp/r0;      % Inner radius of the current distribution (in r0)
        b = (SystemParameters.b)*Rp/r0;   % Outer radius of the current distribution (in r0)
        D = (SystemParameters.D)*Rp/r0;       % Disk Semi-thickness (in r0)
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
        Opacity = 0.2;
        bottom = patch(center(1)+[xout,xin], ...
                   center(2)+[yout,yin], ...
                   z1*ones(1,2*length(xout)),'', 'FaceAlpha', Opacity, 'LineStyle', 'none');
        top = patch(center(1)+[xout,xin], ...
                center(2)+[yout,yin], ...
                z2*ones(1,2*length(xout)),'', 'FaceAlpha', Opacity, 'LineStyle', 'none');
        [X,Y,Z] = cylinder(1,length(xin));
        outer = surf(rout*X+center(1), ...
                 rout*Y+center(2), ...
                 Z*(z2-z1)+z1, 'FaceColor', 'blue', 'FaceAlpha', 0.1, 'LineStyle', 'none');
        inner = surf(rin*X+center(1), ...
                 rin*Y+center(2), ...
                 Z*(z2-z1)+z1, 'FaceColor', 'blue', 'FaceAlpha', 0.1, 'LineStyle', 'none');

        PhiSpanConcatenated = (0:DeltaPhiDeg:PhiMaxDeg)*pi/180;
        ThetaSpanConcatenated = (0:DeltaThetaDeg:ThetaMaxDeg)*pi/180;

            ImgRGB = imread('ColourMap.png');
            [ImRows, ImCols, ImPlanes] = size(ImgRGB);
            ConstructionSurfaceWarp(SurfaceCorrected(:,Jump:Jump:end).', ThetaSpanConcatenated(Jump:Jump:end), PhiSpanConcatenated(1:1:end), ImgRGB)
            set(gcf,'color','w');
            view([45+3*90+5 30])
            axis([-2.3, 2.3, -2.3, 2.3, -2.3, 2.3])
            box on
            grid on

        caxis([-4 0])
        cb = colorbar('eastoutside', 'YTickLabel', {'0.01', '', '0.1', '', '1', '', '10', '', '100'}, 'LineWidth', 1, 'FontName' , 'Helvetica', 'FontSize', 15);
        
        % Optional colorbar positioning for publications
%         cbPosition = cb.Position;
%         ZStretch = 0.8;
%         set(cb, 'Position', [cbPosition(1)*0.99 cbPosition(2)+(cbPosition(4)*(1-ZStretch))/2 cbPosition(3) cbPosition(4)*ZStretch])
            
	hold off
    
end


end