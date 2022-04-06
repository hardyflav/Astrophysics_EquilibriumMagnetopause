

function ParaSystem2 = NoseDetermination(ParaSystem)

% ------------------------------------------------------------------------------------
% Finding nose of the magnetopause, defined as being the closest point of the MP to the Sun
% ------------------------------------------------------------------------------------

ParaSystem2 = ParaSystem;

%% Span Range of Directions and Solve PB in each of them

    % X = [r, theta] relative to Solar direction in NMM
%     NoseFunction = @(X) NosePB(X(1), X(2), ParaSystem);
    NoseFunction = @(X) PB_Nose(X(1), X(2), ParaSystem);

    TiltPhi_Deg = ParaSystem.Tilt.Phi * 180/pi;
    TiltTheta_Deg = ParaSystem.Tilt.Theta * 180/pi;
    Cond_SouthernNose = (0 <= TiltPhi_Deg & TiltPhi_Deg <= 90 | -90 <= TiltPhi_Deg & TiltPhi_Deg <= 0);
    if abs(TiltPhi_Deg) == 90 || TiltTheta_Deg == 0
        AlphaList = 0;
    elseif Cond_SouthernNose
        AlphaList = (-30:0.1:5).' * pi/180;
    else
        AlphaList = (-5:0.1:30).' * pi/180;
    end
        
%     AlphaList = (-30:0.1:30).' * pi/180;
%     AlphaList = (0:0.1:30).' * pi/180;
    rList = 0 * AlphaList;
    
   rGuess = ParaSystem2.rMP .* ParaSystem2.Rp / ParaSystem2.r0;
    for kAlpha = 1:length(AlphaList)
        
        Alpha = AlphaList(kAlpha);
        PBFunctionR = @(r) NoseFunction([r, Alpha]);
        rSol = fsolve(PBFunctionR, rGuess);
        rList(kAlpha) = rSol;
        
    end

%% Finding closest position to the Sun
    if abs(TiltPhi_Deg) == 90 || TiltTheta_Deg == 0
        rMax_r0 = rList;
        AlphaMaxR_Rad = AlphaList;
        AlphaMax_Rad = AlphaList;
        XMax_r0 = rMax_r0;
    else
        InterpMethod = 'linear';
        rInterpolant = griddedInterpolant(AlphaList, rList, InterpMethod);
        Xinterpolant = griddedInterpolant(AlphaList, rList .* cos(AlphaList), InterpMethod);
        AlphaList_Eval = (AlphaList(1)*180/pi :0.1:AlphaList(end)*180/pi).'  * pi/180;
        rList_Eval = rInterpolant(AlphaList_Eval);
        XList_Eval = Xinterpolant(AlphaList_Eval);


        rMax_r0 = max(rList_Eval);
        AlphaMaxR_Rad = mean(AlphaList_Eval( rList_Eval == rMax_r0 ));

        XMax_r0 = max(XList_Eval);
        AlphaMax_Rad = mean(AlphaList_Eval( XList_Eval == XMax_r0 ));

        NoseFunction([rMax_r0; AlphaMaxR_Rad]);
        NoseFunction([XMax_r0 / cos(AlphaMax_Rad); AlphaMax_Rad]);
    end
    %{
    % Figure: X(alpha) = f(alpha), to explain how the position of the Nose
    % is determined
    figure;
    hold on
        Opacity = 0.5;
        PX = plot( AlphaList_Eval .* 180/pi,  XList_Eval .* ParaSystem.r0 ./ ParaSystem.Rp );
        PMax_X = scatter(AlphaMax_Rad .* 180/pi, XMax_r0 .* ParaSystem.r0 ./ ParaSystem.Rp , 40, 'filled');
        PLine = line([AlphaMax_Rad .* 180/pi AlphaMax_Rad .* 180/pi], [0 XMax_r0 .* ParaSystem.r0 ./ ParaSystem.Rp]);
        Colour_PMax_X = cbrewer2('Purples', 20);
        set([PX], 'LineStyle', '-','linewidth', 6, 'Color', horzcat( Colour_PMax_X(17, :), Opacity));
        set(PMax_X, 'MarkerEdgeColor', 'black', 'LineWidth', 1, 'SizeData', 200, 'CData', Colour_PMax_X(17, :));
        Ylabel1 = ylabel('X_0 / R_S (Full lines)');
        axes = gca;
        set(axes, 'FontName'   , 'Helvetica',...
            'XAxisLocation', 'bottom',...
            'FontSize', 16,...
            'LineWidth'   , 1 ...
            );
        axis([-30, 0, 0, 20])
    
        yyaxis right
        Pr = plot( AlphaList_Eval .* 180/pi,  rList_Eval .* ParaSystem.r0 ./ ParaSystem.Rp );
        PMax_r = scatter(AlphaMaxR_Rad .* 180/pi, rMax_r0 .* ParaSystem.r0 ./ ParaSystem.Rp , 40, 'filled');
%         Colour_PMax_X = cbrewer2('Greys', 20);
        set([Pr], 'LineStyle', '-.', 'linewidth', 6, 'Color', horzcat( Colour_PMax_X(13, :), Opacity) );
        set([PX], 'LineStyle', '-', 'linewidth', 6, 'Color', horzcat( Colour_PMax_X(13, :), Opacity) );
        set(PMax_r,  'MarkerEdgeColor', 'black', 'LineWidth', 1, 'SizeData', 200, 'CData', Colour_PMax_X(13, :));
        set(PLine, 'linestyle', '-.', 'linewidth', 2, 'Color', Colour_PMax_X(17, :))
    
    
    % ---- Plot options
        axes = gca;
        set(axes, 'FontName'   , 'Helvetica',...
            'XAxisLocation', 'bottom',...
            'FontSize', 16,...
            'LineWidth'   , 1 ...
            );
        Ylabel = ylabel('r_0 / R_S (Dashed lines)');
        Xlabel = xlabel('\alpha (\circ)');
        set([Xlabel, Ylabel1, Ylabel], 'FontName', 'AvantGarde', 'FontSize', 18);
        axes.YColor = 'black' 
        grid on
        box on
        set(gcf,'color','w');
        ColourBackground = [0, 0, 0, 0.3];
        set(gca,'Color',ColourBackground)
    
        axis([-20, 20, 20, 26])
    
    
    export_fig DeterminationPositionNose.pdf -opengl -q101
    %}
    
    
    
%% Computing corresponding coordinates of the nose

%     rMax_Rp = rMax_r0 .* ParaSystem.r0 ./ ParaSystem.Rp;
    rMax_Rp = (XMax_r0 .* ParaSystem.r0 ./ ParaSystem.Rp) / cos(AlphaMax_Rad);
    AlphaMax_Deg = AlphaMax_Rad * 180/pi;
    
%     hold on
%     scatter(ParaSystem.Tilt.Phi * 180/pi, AlphaMax_Deg, 'filled');
    
    if (AlphaMax_Rad)>0
        NosePhi = (pi/2-ParaSystem.AlphaDipole);
        AlphaMax_Rad = abs(AlphaMax_Rad);
    elseif (AlphaMax_Rad)<0
        NosePhi = -(pi/2+ParaSystem.AlphaDipole);
        AlphaMax_Rad = abs(AlphaMax_Rad);
    else
        NosePhi = (pi/2-ParaSystem.AlphaDipole);
    end

    X_KSMU = rMax_Rp*cos(AlphaMax_Rad);
    Y_KSMU = rMax_Rp*sin(AlphaMax_Rad)*cos(NosePhi);
    Z_KSMU = rMax_Rp*sin(AlphaMax_Rad)*sin(NosePhi);
    
    NosePosition_KSMU = [X_KSMU; Y_KSMU; Z_KSMU];
    Nose.KSMU_Rp = NosePosition_KSMU;

%% Plotting Nose and its Direction 
    Radius = 0.4;
    [u,v,w]=sphere(25);
    NosePlot = surf(NosePosition_KSMU(1)+Radius*u,NosePosition_KSMU(2)+Radius*v,NosePosition_KSMU(3)+Radius*w);
    set(NosePlot, 'FaceLighting', 'none', 'FaceColor', 'red', 'EdgeColor', 'none', 'FaceAlpha', 0.5)
    
    NoseDirection = line( [0 NosePosition_KSMU(1)], [0, NosePosition_KSMU(2)], [0, NosePosition_KSMU(3)] );
    Colour_Equator = cbrewer2('Greens', 10);
    set(NoseDirection, 'Color', horzcat(Colour_Equator(8,:), 0.3), 'linewidth', 2, 'LineStyle', '-.')


ParaSystem2.Nose = Nose;
Nose_r0 = ParaSystem2.Nose.KSMU_Rp .* ParaSystem.Rp ./ ParaSystem.r0;
ParaSystem2.NoseDirection = Nose_r0 ./ norm(Nose_r0);


%% Computing (Nose -> KSM) Conversion Matrix
    exNose = ParaSystem2.NoseDirection;
    
    MTilted = ParaSystem2.Tilt.RotationPhi*( ParaSystem2.Tilt.RotationTheta* [0; 0; -1]);
    ezNose_temp = -(MTilted - dot(exNose, MTilted)*exNose);
    ezNose = ezNose_temp ./ norm(ezNose_temp);
    
    eyNose = cross(ezNose, exNose);
    
    ParaSystem2.P_CartNose2KSM  = [exNose, eyNose, ezNose];
    
end
