


AU = 149597871000;
t = linspace(0, 10^7, 1000 );
u = 400 * 1000 ;
r0 = 0;

r = r0 + u*t;
omega = 2.7 * 10^(-6);
phi0 = pi;
phi = -omega*t + phi0;
x = r .* cos(phi) ./ AU;
y = r .* sin(phi) ./ AU;
figure;
hold on
    Spiral = plot(x, y);
    Col_Spiral = cbrewer2('Oranges', 20);
    set(Spiral, 'LineWidth', 7, 'Color', [Col_Spiral(18,:), 0.4])
    axis equal
    axis([-12 12 -12 12])
    
    rMer = 0.387;
    rEarth = 1;
    rJupiter = 5.2;
    rSaturn = 9.5;
    Col_Orbits = cbrewer2('BuGn', 20);
    
    
    rList = [rMer, rEarth, rJupiter, rSaturn];
    theta = linspace(0, 2*pi, 100);
    Plot_Mer = plot( rList(1) .* cos(theta), rList(1) .* sin(theta));
    Plot_Earth = plot( rList(2) .* cos(theta), rList(2) .* sin(theta));
    Plot_Jup = plot( rList(3) .* cos(theta), rList(3) .* sin(theta));
    Plot_Sat = plot( rList(4) .* cos(theta), rList(4) .* sin(theta));
    
    set(Plot_Mer, 'LineWidth', 3, 'LineStyle', '-.', 'Color', [Col_Orbits(13,:), 0.5]);
    set(Plot_Earth, 'LineWidth', 3, 'LineStyle', '-.', 'Color', [Col_Orbits(15,:), 0.5]);
    set(Plot_Jup, 'LineWidth', 3, 'LineStyle', '-.', 'Color', [Col_Orbits(17,:), 0.5])
    set(Plot_Sat, 'LineWidth', 3, 'LineStyle', '-.', 'Color', [Col_Orbits(19,:), 0.5]);
    
    axes = gca;
    set(axes, 'FontName'   , 'Helvetica',...
        'XAxisLocation', 'bottom',...
        'YAxisLocation', 'left',...
        'FontSize', 16,...
        'LineWidth'   , 1 ...
        );
    Ylabel = ylabel('Distance from Sun (AU)');
    Xlabel = xlabel('Distance from Sun (AU)');
    ColourBackground = [0, 0, 0, 0.3];
    set(gca,'Color',ColourBackground)
    set(gcf,'color','w');
    grid on
    
    export_fig ParkerSpiral.pdf -opengl -q101
    