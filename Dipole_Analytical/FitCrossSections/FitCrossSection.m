


% ThetaVector = horzcat(5:15:65, 70:5:80, 100:20:120);
% ListCoefficients = ones(length(ThetaVector), 5);
% 
% C = 1;
% for i = ThetaVector
%     ListCoefficients(C,:) = FitCrossSection(SurfaceCorrected, 1/2, i);
%     C = C+1;
% end
% 
% fliplr(ListCoefficients)
% 
% ThetaVector = horzcat(1:5:120);
% ListCoefficients = ones(length(ThetaVector), 5);
% 
% C = 1;
% for i = ThetaVector
%     ListCoefficients(C,:) = FitCrossSection(SurfaceCorrected, 1/2, i);
%     C = C+1;
% end
% 
% fliplr(ListCoefficients)
% 


function [PolynomialFit] = FitCrossSection(Surface, DeltaThetaDeg, ThetaValueDeg)

ThetaPositionGrid = ThetaValueDeg/DeltaThetaDeg + 1;
rValuesCrossSection = Surface(:, ThetaPositionGrid);

RhoValues = rValuesCrossSection.*sin(ThetaValueDeg*pi/180);
PhiValues = (0:DeltaThetaDeg:90)*pi/180;

PolynomialFit = polyfit(PhiValues(1:2:end), RhoValues(1:2:end).', 4);
FunctionFit = polyval(PolynomialFit,PhiValues);

hold on
[X, Y] = pol2cart(PhiValues, RhoValues.');
plot(X,Y, '-', 'linewidth', 4)



[XFit, YFit] = pol2cart(PhiValues, FunctionFit);
plot(XFit(1:10:end), YFit(1:10:end), '.', 'MarkerSize', 16)

ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
axis([0 1.7 0 1.7])
axis square
grid on
end

