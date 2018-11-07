function ConstructionSurfaceWarp = ConstructionSurfaceWarp(Array, theta_span, phi_span, ImgRGB)


angles_x = arrayfun(@sin, theta_span(1:end) ).'*arrayfun(@cos,  phi_span(1:end) );
angles_y = arrayfun(@sin, theta_span(1:end) ).'*arrayfun(@sin, phi_span(1:end) );
angles_z = arrayfun(@cos, theta_span(1:end) ).'*ones(1,size(phi_span,2));

X = Array .* angles_x;
Y = Array .* angles_y;
Z = Array .* angles_z;

% Surface
% hold on
%      C = Z;
%      surf(X,Y,Z, C,'FaceAlpha', 0.5);
%      surf(X,-Y,Z, C,'FaceAlpha', 0.5);
%      surf(-X,Y,Z, C,'FaceAlpha', 0.5);
%      surf(-X,-Y,Z, C,'FaceAlpha', 0.5);
%       
%     xlabel('X','FontSize',14, 'FontWeight','bold');
%     ylabel('Y','FontSize',14, 'FontWeight','bold');
%     zlabel('Z','FontSize',14, 'FontWeight','bold');
%     set(gca, 'FontSize', 14)
%     axis square
%     grid on  
%     box on
%     
%     mArrow3([0,1/2,0],[0,-1/2,0],'color','red');
%     
% hold off

% Texture
hold on
     
     Texture1 = warp(X, Y, Z, ImgRGB);
     Texture2 = warp(X, -Y, Z, ImgRGB);  
     Texture3 = warp(-X, Y, Z, ImgRGB);  
     Texture4 = warp(-X, -Y, Z, ImgRGB);  
     
     Opacity = 0.7*1;
     set(Texture1, 'FaceAlpha', Opacity);
     set(Texture2, 'FaceAlpha', Opacity);
     set(Texture3, 'FaceAlpha', Opacity);
     set(Texture4, 'FaceAlpha', Opacity);
     
%     xlabel('X','FontSize',14, 'FontWeight','bold');
%     ylabel('Y','FontSize',14, 'FontWeight','bold');
%     zlabel('Z','FontSize',14, 'FontWeight','bold');
%     set(gca, 'FontSize', 14)
%     axis equal
%     grid on  
%     box on
%     set(gcf,'color','w');
    
    
        axes = gca;
        set(axes, 'FontName'   , 'Helvetica',...
                  'XAxisLocation', 'origin',...
                  'YAxisLocation', 'origin',...
                  'FontSize', 16,...
                  'LineWidth'   , 1 ...
              );
        Ylabel = ylabel('X (r_0)');
        Xlabel = xlabel('Z (r_0)');  
        Zlabel = zlabel('Y (r_0)');  
        
        set([Xlabel, Ylabel, Zlabel], 'FontName'   , 'AvantGarde', 'FontSize', 18);
%         set(Ylabel, 'position', [-0.45 1.35 0]);
        axis equal
        grid on
%         box on
        set(gcf,'color','w')      
    
    
    
    view([0 90])    
    
    
    
    
    
    mArrow3([0,-1/3,0],[0,1/3,0],'color','red','stemWidth', 0.01);
    
              
    Bp = 20000*10^(-9);     % Equatorial field, T
    Rp = 60280*10^3;        % Planet radius, m
    mu_0 = 4*pi*10^(-7);
    Pswnpa = 0.02;
    b0 = sqrt(2*mu_0*Pswnpa*10^(-9));   % Field scale
    r0 = (2*Bp*Rp^3/b0)^(1/3);          % Distance scale

    RpScaled = Rp/r0;
    [Xs, Ys, Zs] = sphere;
    surf(RpScaled*Xs, RpScaled*Ys, RpScaled*Zs)
    
    
hold off




% 
% hold on
% 
%    % Surface
%       C = Z;
%       surf(X,Y,Z, C,'FaceAlpha', 0.);
%       surf(X,-Y,Z, C,'FaceAlpha', 0.);
%       surf(-X,Y,Z, C,'FaceAlpha', 0.);
%       surf(-X,-Y,Z, C,'FaceAlpha', 0.);
%       
%    % Texture
%    
%      Texture1 = warp(-X, Y, Z, ImgRGB);  
%      Texture2 = warp(X, -Y, Z, ImgRGB);  
%      Texture3 = warp(-X, -Y, Z, ImgRGB);  
%      Texture4 = warp(X, Y, Z, ImgRGB);
%      
%      set(Texture1, 'FaceAlpha', 0.9);
%      set(Texture2, 'FaceAlpha', 0.9);
%      set(Texture3, 'FaceAlpha', 0.9);
%      set(Texture4, 'FaceAlpha', 0.9);
%      
%    % Lighting
% %     light('Position',[0 0 0.05],'Style','infinite');
% 
%    %Axes
%     xlabel('X','FontSize',14,'FontWeight','bold');
%     ylabel('Y','FontSize',14,'FontWeight','bold');
%     zlabel('Z','FontSize',14,'FontWeight','bold');
%     set(gca, 'FontSize', 14)
%     axis square
%     grid on  
%     box on
%       
%     mArrow3([0,1/2,0],[0,-1/2,0],'color','red');
%      
%       
% hold off




% hold on
%     contour3(X, Y, Z, 40, '.', 'linewidth', 2);
%     contour3(X, -Y, Z, 40, '.', 'linewidth', 2);
%     contour3(-X, Y, Z, 40, '.', 'linewidth', 2);
%     contour3(-X, -Y, Z, 40, '.', 'linewidth', 2);
%     axis equal
%     grid on
%     box on
% hold off

end

