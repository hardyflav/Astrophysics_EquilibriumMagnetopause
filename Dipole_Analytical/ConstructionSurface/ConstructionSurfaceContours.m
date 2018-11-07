function ConstructionSurfaceContours = ConstructionSurfaceContours(Array, theta_span, phi_span)


angles_x = arrayfun(@sin, theta_span(1:end) ).'*arrayfun(@cos,  phi_span(1:end) );
angles_y = arrayfun(@sin, theta_span(1:end) ).'*arrayfun(@sin, phi_span(1:end) );
angles_z = arrayfun(@cos, theta_span(1:end) ).'*ones(1,size(phi_span,2));

X = Array .* angles_x;
Y = Array .* angles_y;
Z = Array .* angles_z;


    hold on

      C = Z;
%       surf(X,Y,Z, C,'FaceAlpha', 0.3);
%       surf(X,-Y,Z, C,'FaceAlpha', 0.3);
%       surf(-X,Y,Z, C,'FaceAlpha', 0.3);
%       surf(-X,-Y,Z, C,'FaceAlpha', 0.3);

%       surf(Z, X,Y, C,'FaceAlpha', 0.3);
%       surf(Z, X,-Y, C,'FaceAlpha', 0.3);
%       surf(Z, -X,Y, C,'FaceAlpha', 0.3);
%       surf(Z, -X,-Y, C,'FaceAlpha', 0.3);


      xlabel('X','FontSize',14,'FontWeight','bold');
      ylabel('Y','FontSize',14,'FontWeight','bold');
      zlabel('Z','FontSize',14,'FontWeight','bold');
      set(gca, 'FontSize', 14)
        axis equal
%         shading interp
      
%       axis([-1,1, -0.9,0.9, -1.5,1.3])


%       mArrow3([0, 0,1/2],[0, 0,-1/2],'color','red');
      mArrow3([0,1/3,0],[0,-1/3,0],'color','red','stemWidth', 0.01);
      
    contourf(X, Y, Z, 20, '-k', 'linewidth', 1);
    contourf(X, -Y, Z, 20, '-k', 'linewidth', 1);
    contourf(-X, Y, Z, 20, '-k', 'linewidth', 1);
    contourf(-X, -Y, Z, 20, '-k', 'linewidth', 1);
    
      box on
      
          
    Bp = 20000*10^(-9);     % Equatorial field, T
    Rp = 60280*10^3;        % Planet radius, m
    mu_0 = 4*pi*10^(-7);
    Pswnpa = 0.02;
    b0 = sqrt(2*mu_0*Pswnpa*10^(-9));   % Field scale
    r0 = (2*Bp*Rp^3/b0)^(1/3);          % Distance scale

    RpScaled = Rp/r0;
    [Xs, Ys, Zs] = sphere;
    surf(RpScaled*Xs, RpScaled*Ys, RpScaled*Zs)
    
    set(gcf,'color','w');
    grid on
    hold off

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

