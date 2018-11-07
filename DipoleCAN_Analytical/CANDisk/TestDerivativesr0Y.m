


% Offset = (-1:0.001:1);
% theta0Y_Before = (90-Offset).*pi/180;
% theta0Y = pi/2;
% theta0Y_After = (90+Offset).*pi/180;
rVal = r0Y;
% 
% drMeridianCAN(theta0Y_Before, rVal, beta, CAN_DiskParameters)
% drMeridianCAN(theta0Y, rVal, beta, CAN_DiskParameters)
% drMeridianCAN(theta0Y_After, rVal, beta, CAN_DiskParameters)
% 
% 
% drMeridianCAN2(theta0Y_Before, rVal, beta, CAN_DiskParameters)
% drMeridianCAN2(theta0Y, rVal, beta, CAN_DiskParameters)
% drMeridianCAN2(theta0Y_After, rVal, beta, CAN_DiskParameters)


Offset = (-1:0.01:1);
theta0Y_After = (90+Offset).*pi/180;
X = Offset;
Y1=[];
Y2=[];

for k = 1:length(Offset)
    Y1 = horzcat(Y1, drMeridianCAN(theta0Y_After(k), rVal, beta, CAN_DiskParameters));
	Y2 = horzcat(Y2, drMeridianCAN2(theta0Y_After(k), rVal, beta, CAN_DiskParameters));
end

hold on
plot(X, Y1, 'o', 'LineWidth', 3)
plot(X, Y2, 'o', 'LineWidth', 3)
axis([-1, 1, -1, 0.6])

