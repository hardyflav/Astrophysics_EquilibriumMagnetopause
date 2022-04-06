
%% Import Data
DegToRad = pi/180;
PathFile = '/Users/flavienhardy/Documents/Git/MagnetopauseTilt/Plots/CassiniMission/data/cassini_orbit.txt';
File = importdata(PathFile);
Time.Year = File.data(:, 1);
Time.Day = File.data(:, 2);
Orbit.X = File.data(:, 6);
Orbit.Y = File.data(:, 7);
Orbit.Z = File.data(:, 8);

Time.Proxy = (Time.Year+Time.Day/1000);
Orbit.XInterp = griddedInterpolant( Time.Proxy,  Orbit.X, 'spline' );
Orbit.YInterp = griddedInterpolant( Time.Proxy,  Orbit.Y, 'spline' );
Orbit.ZInterp = griddedInterpolant( Time.Proxy,  Orbit.Z, 'spline' );


%% Plot Orbit
EquinoxMission = (Time.Proxy >= 2004 & Time.Year <= 2006);
EquinoxMission_Time_temp = Time.Proxy(EquinoxMission);
EquinoxMission_Time = ( EquinoxMission_Time_temp(1) : 1/3600 : EquinoxMission_Time_temp(end) );

X = Orbit.XInterp(EquinoxMission_Time);
Y = Orbit.YInterp(EquinoxMission_Time);
Z = Orbit.ZInterp(EquinoxMission_Time);
figure;
hold on
    PlotOrbit = plot3( X(1:1:end), -Y(1:1:end), Z(1:1:end) );
    ColorOrbit = cbrewer2('Pastel1', 9);
    Opacity = 0.3;
    set(PlotOrbit, 'linewidth', 2, 'color', horzcat(ColorOrbit(3, :), Opacity) );
    axis equal
    
copyobj(gcf, 0)
figure
     PlotOrbit = plot3(Orbit.X(EquinoxMission), Orbit.Y(EquinoxMission), Orbit.Z(EquinoxMission))
    set(PlotOrbit, 'linewidth', 2, 'color', horzcat(ColorOrbit(3, :), Opacity) );

%% Plot Surface
figure;
        hold on
        ColorSurface = cbrewer2('PuBu', 20);
        Surface = surf(X_KSM, Y_KSM, Z_KSM, X_KSM);

        YNose = -YNose;           
        X_KSM = P_CartNose2KSM(1,1) .* XNose + P_CartNose2KSM(1,2) .* YNose + P_CartNose2KSM(1,3) .* ZNose;
        Y_KSM = P_CartNose2KSM(2,1) .* XNose + P_CartNose2KSM(2,2) .* YNose + P_CartNose2KSM(2,3) .* ZNose;
        Z_KSM = P_CartNose2KSM(3,1) .* XNose + P_CartNose2KSM(3,2) .* YNose + P_CartNose2KSM(3,3) .* ZNose;
        Surface2 = surf(X_KSM, Y_KSM, Z_KSM, X_KSM);        
        
        shading flat
        set([Surface, Surface2], 'FaceAlpha', 0.6, 'EdgeColor', 'none');
        axis equal
        axis([-50 50 -25 25 -25 25])
                
        colormap(cbrewer2('BuPu', 30))

        
%% Axes And BackGround
    ScaleRp = 30;
    Opacity = 0.5;
    ColorLine = [1, 1, 1, Opacity];
    % ---- Axes at Origin
        L1 = plot3([-ScaleRp ScaleRp],[0 0],[0 0],'Color', ColorLine,...
            'LineWidth', 2,...
            'LineStyle', "-.");

        L2 = plot3([0 0],[-ScaleRp ScaleRp],[0 0],'Color', ColorLine,...
            'LineWidth', 2,...
            'LineStyle', "-.");
        
        L3 = plot3([0 0],[0 0],[-ScaleRp ScaleRp],'Color', ColorLine,...
            'LineWidth', 2,...
            'LineStyle', "-.");
        
        ColourBackground = [0, 0, 0, 0.5];
        set(gca,'Color',ColourBackground)
       
        box off
        grid off
        
        Ylabel = ylabel('X (R_p)');
        Xlabel = xlabel('Z (R_p)');
        Zlabel = zlabel('Y (R_p)');
        
        axis([-25 35 -35 35 -35 35])

        
        
        
        
        
        
        
        
%% ON THE SIDE

%{

        Cond = (abs(YNose)<0.1 & ZNose > 0);
        NMM = plot3(X_KSM(Cond), Y_KSM(Cond), Z_KSM(Cond));
        set(NMM, 'LineWidth', 2)
