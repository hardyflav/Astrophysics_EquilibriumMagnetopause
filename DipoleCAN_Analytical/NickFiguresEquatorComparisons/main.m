

figure(1);    

%% List PLots
ThetaAxis = ThetaEquator;

[YDipole, ZDipole] = pol2cart(ThetaAxis, rEquator_Dipole);

[YDipoleDisk, ZDipoleDisk] = pol2cart(ThetaAxis, rEquator_DipoleDisk);

[YDipoleDiskBeta10, ZDipoleDiskBeta10] = pol2cart(ThetaAxis, rEquator_DipoleDiskBeta10);
[YDipoleDiskBeta25, ZDipoleDiskBeta25] = pol2cart(ThetaAxis, rEquator_DipoleDiskBeta25);

X = zeros(length(YDipole), 4);
X(:,1) = YDipole;
X(:,2) = YDipoleDisk;
X(:,3) = YDipoleDiskBeta10;
X(:,4) = YDipoleDiskBeta25;

Y = zeros(length(ZDipole), 4);
Y(:,1) = ZDipole;
Y(:,2) = ZDipoleDisk;
Y(:,3) = ZDipoleDiskBeta10;
Y(:,4) = ZDipoleDiskBeta25;


Xline = [1, 1];
Yline = [0, 3];

X1 = -1;
X2 = 2;
Y1 = 0;
Y2 = 2.5;



%% Plotting


[ha, ~] = tight_subplot(3,1,[.01 .03],[.1 .1],[.05 .05]);
set(ha(1:3), 'YTickLabelMode', 'auto',...
             'Xdir', 'reverse',...
             'YTickLabel', {'0', '', '1', '', '2', ''},...
             'YAxisLocation', 'left',...
             'FontSize', 16,...
             'LineWidth'   , 1);
set(ha(3), 'XTickLabelMode', 'auto') 
set(ha(1), 'XAxisLocation', 'top', 'XTickLabelMode', 'auto') 


% Subfigure 1: Dipole
axes(ha(1));
Ylabel = ylabel('Y (r_0)');
Xlabel = xlabel('X (r_0)');
set( gca                       , ...
    'FontName'   , 'Helvetica');
set([Xlabel, Ylabel], ...
    'FontName'   , 'AvantGarde');
axis equal
grid on
text(1.8,2.3,'a)',...
    'FontSize',16,'HorizontalAlignment','left')
hold on
        Dipole = plot(X(:,1), Y(:,1), 'LineWidth', 4);
        plot(Xline, Yline, '--k', 'LineWidth', 1) 
        axis([X1 X2 Y1 Y2]) 
hold off
set(Dipole, 'Color', [0.8500    0.3300    0.1000]);



% Subfigure 2: Dipole + Disk
axes(ha(2))
Ylabel = ylabel('Y (r_0)');
set( gca                       , ...
    'FontName'   , 'Helvetica' );
set([Ylabel], ...
    'FontName'   , 'AvantGarde');axis equal
grid on
text(1.8,2.3,'b)',...
    'FontSize',16,'HorizontalAlignment','left')
hold on;
        DipoleDisk = plot(X(:,2), Y(:,2), 'LineWidth', 4);
        plot(Xline, Yline, '--k', 'LineWidth', 1) 
        axis([X1 X2 Y1 Y2])  
hold off
set(DipoleDisk, 'Color', [0.4700    0.6700    0.1900]);



% Subfigure 3: Dipole + Disk + Beta
axes(ha(3));
axis equal
Ylabel = ylabel('Y (r_0)');
Xlabel = xlabel('X (r_0)');
set( gca                       , ...
    'FontName'   , 'Helvetica' );
set([Xlabel, Ylabel], ...
    'FontName'   , 'AvantGarde');
grid on
text(1.8,2.3,'c)',...
    'FontSize',16,'HorizontalAlignment','left')
hold on
        beta10 = plot(X(:,3), Y(:,3), 'LineWidth', 4);
        beta25 = plot(X(:,4), Y(:,4), 'LineWidth', 4);
        plot(Xline, Yline, '--k', 'LineWidth', 1) 
        axis([X1 X2 Y1 Y2])    
hold off

Legend = legend([beta10, beta25],...
    '\beta = 10',...
    '\beta = 25',...
    'location', 'SouthEast');
set([Legend], ...
    'FontSize'   , 16           );

set(beta10, 'Color', [.75 .75 1]);
set(beta25, 'Color', [0    0.4500    0.7400]);




set(gcf,'renderer','Painters')
print -depsc2 -tiff -r300 -painters EquatorialSolutions.eps





% c = uisetcolor([0.6 0.8 1])




