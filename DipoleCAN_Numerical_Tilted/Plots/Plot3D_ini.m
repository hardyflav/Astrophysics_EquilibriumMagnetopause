
function Plot = Plot3D_ini(Rotation, ParaSystem)

    Plot = figure;
        ScaleRp = 30;

        
    % ---- Dipole Field Lines
        ColourFieldLines = cbrewer2('Reds', 10);
        FieldLines = lforce3d(10, 20, 'k', ParaSystem, Rotation, ColourFieldLines(8, :));
        
    % ---- Planet and Rings
        Radius = 1;
        [u,v,w]=sphere(25);
        Sphere = surf(Radius*u,Radius*v,Radius*w);
        colormap('default');
        camlight right;
        lighting phong;
        ColorPlanet = cbrewer2('YlOrBr', 10);
        set(Sphere, 'EdgeColor', 'none', 'FaceColor', ColorPlanet(4,:))
        
        r = 66900*10^3/ParaSystem.Rp; %inner radius
        R = 140180*10^3/ParaSystem.Rp; %outer radius
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
        
    % shading interp

        if strcmp(Rotation, 'DipoleTilted')
                rotate( [Sphere, Rings] , [0, 1, 0], ParaSystem.Tilt.Theta  * 180/pi, [0, 0, 0])
                rotate( [Sphere, Rings] , [0, 0, 1], ParaSystem.Tilt.Phi * 180/pi, [0, 0, 0])
        end
        
        
        
    % ---- Plane Plot options
        Colour_NMM = cbrewer2('Oranges', 10);
        Colour_Equator = cbrewer2('Greens', 10);

    % ---- Solar Wind
        ArrowSolarWind = mArrow3([25; 0; 10], [15; 0; 10],'color','red','stemWidth', 0.3);
        set (ArrowSolarWind, 'FaceAlpha', 0.4)
        if strcmp(Rotation, 'DipoleAligned')
            rotate( ArrowSolarWind , [0, 0, 1], -ParaSystem.Tilt.Phi * 180/pi, [0, 0, 0]);
            rotate( ArrowSolarWind , [0, 1, 0], -ParaSystem.Tilt.Theta  * 180/pi, [0, 0, 0]);
        end
        
    % ---- Magnetic Dipole
        ArrowDipole = mArrow3([0; 0; 7], [0; 0; -7],'color','red','stemWidth', 0.2);
        if strcmp(Rotation, 'DipoleTilted')
            rotate( ArrowDipole , [0, 1, 0], ParaSystem.Tilt.Theta  * 180/pi, [0, 0, 0]);
            rotate( ArrowDipole , [0, 0, 1], ParaSystem.Tilt.Phi * 180/pi, [0, 0, 0]);
        end
        
    % ---- NMM Plane
        AlphaDipole = ParaSystem.AlphaDipole;
        
        XList = [-80 80 80 -80];
        ZList = [-80 -80 80 80];
        YList = [0 0 0 0];
            
        NMM_Plane = patch('XData',XList,'YData',YList,'ZData',ZList);
        rotate( NMM_Plane , [1, 0, 0], -AlphaDipole * 180/pi, [0, 0, 0]);
        
        Color_NMMPlane = cbrewer2('Greens', 40);
        Opacity = 0.5;
        set(NMM_Plane, 'FaceColor', Color_NMMPlane(10,:), 'FaceAlpha', Opacity, 'EdgeColor', 'none');
        
        if strcmp(Rotation, 'DipoleAligned')
            rotate( NMM_Plane , [0, 0, 1], ParaSystem.Tilt.Phi * 180/pi, [0, 0, 0]);
            rotate( NMM_Plane , [0, 1, 0], ParaSystem.Tilt.Theta  * 180/pi, [0, 0, 0]);
        end
        

    % ---- Magnetic Equator Plane

        XList = [-80 80 80 -80];
        YList = [-80 -80 80 80];
        ZList = [0 0 0 0];
        
        
        Equ_Plane = patch('XData',XList,'YData',YList,'ZData',ZList);
        rotate( Equ_Plane , [0, 1, 0], ParaSystem.Tilt.Theta  * 180/pi, [0, 0, 0]);
        rotate( Equ_Plane , [0, 0, 1], ParaSystem.Tilt.Phi  * 180/pi, [0, 0, 0]);
        
        Color_EquPlane = cbrewer2('Oranges', 40);
        Opacity = 0.5;
        set(Equ_Plane, 'FaceColor', Color_EquPlane(10,:), 'FaceAlpha', Opacity, 'EdgeColor', 'none');
        
        if strcmp(Rotation, 'DipoleAligned')
            rotate( Equ_Plane , [0, 0, 1], ParaSystem.Tilt.Phi * 180/pi, [0, 0, 0]);
            rotate( Equ_Plane , [0, 1, 0], ParaSystem.Tilt.Theta  * 180/pi, [0, 0, 0]);
        end
        
        
        % ---- Plot options
        axes = gca;
        set(axes, 'FontName'   , 'Helvetica',...
            'XAxisLocation', 'top',...
            'YAxisLocation', 'left',...
            'FontSize', 16,...
            'LineWidth'   , 1 ...
            );
        Ylabel = ylabel('X (R_p)');
        Xlabel = xlabel('Z (R_p)');
        Zlabel = zlabel('Y (R_p)');

        set([Xlabel, Ylabel, Zlabel], 'FontName', 'AvantGarde', 'FontSize', 18);
        axis equal
        grid on
        box on
        set(gcf,'color','w');
        ColourBackground = [0, 0, 0, 0.3];
        set(gca,'Color',ColourBackground)
        axis([-10, ScaleRp, -ScaleRp, ScaleRp, -ScaleRp*1/2, ScaleRp*3/4])
        view(110, 20)
        
    % ---- Axes at Origin
        plot3([-ScaleRp ScaleRp],[0 0],[0 0],'Color', [0.5020    0.5020    0.5020],...
            'LineWidth', 1,...
            'LineStyle', "-.");

        plot3([0 0],[-ScaleRp ScaleRp],[0 0],'Color', [0.5020    0.5020    0.5020],...
            'LineWidth', 1,...
            'LineStyle', "-.");
        plot3([0 0],[0 0],[-ScaleRp ScaleRp],'Color', [0.5020    0.5020    0.5020],...
            'LineWidth', 1,...
            'LineStyle', "-.");



end