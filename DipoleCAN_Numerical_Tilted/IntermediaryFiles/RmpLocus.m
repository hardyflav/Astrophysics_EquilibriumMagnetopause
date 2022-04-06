
% Psw Fixed at 0.01 nPa

% -------------- Dawn Side --------------------
PhiList_Dawn = horzcat( 0, 10, 20, 30, 40:20:80, 90).' ;

NoseValues_Dawn = [
                    [ 24.7428    0.0000   -6.8153 ];...
                    [ 24.7764   -0.5270   -6.6644 ];...
                    [ 24.9211   -0.9935   -6.3210 ];...
                    [ 25.1261   -1.3588   -5.8296 ];...
                    [ 25.3040   -1.5696   -5.1444 ];...
                    [ 25.6857   -1.3654   -3.1932 ];...
                    [ 26.0567   -0.5494   -1.0991 ];...
                    [ 26.0988         0         0 ];...
                  ];

 % Building Phi=(0:180) From Phi=(0:90)
 PhiList_Dawn_Tot = horzcat( 0:10:40, 60:20:80, 90, 100:20:120, 140, 150, 160, 170, 180).' ;
 NoseValues_Dawn_Winter = flipud(NoseValues_Dawn(1:end-1, :));
 NoseValues_Dawn_Winter(:, 2) =  -NoseValues_Dawn_Winter(:, 2);
 NoseValues_Dawn_Winter(:, 3) =  -NoseValues_Dawn_Winter(:, 3);

 NoseValues_Dawn_Tot = vertcat( NoseValues_Dawn, NoseValues_Dawn_Winter );


              
% -------------- Dusk Side --------------------

PhiList_Dusk_Tot = -PhiList_Dawn_Tot ;
           
NoseValues_Dusk_Tot = NoseValues_Dawn_Tot;
NoseValues_Dusk_Tot(:, 2) = -NoseValues_Dawn_Tot(:, 2);

              
% -------------- Concatenation --------------

PhiList = vertcat( flipud(PhiList_Dusk_Tot), PhiList_Dawn_Tot(2:end) );
NoseValues_X = vertcat( flipud(NoseValues_Dusk_Tot(:,1)), NoseValues_Dawn_Tot(2:end,1) );
NoseValues_Y = vertcat( flipud(NoseValues_Dusk_Tot(:,2)), NoseValues_Dawn_Tot(2:end,2) );
NoseValues_Z = vertcat( flipud(NoseValues_Dusk_Tot(:,3)), NoseValues_Dawn_Tot(2:end,3) );

Method = 'spline';
NoseValues_XInterp = griddedInterpolant( PhiList, NoseValues_X, Method );
NoseValues_YInterp = griddedInterpolant( PhiList, NoseValues_Y, Method );
NoseValues_ZInterp = griddedInterpolant( PhiList, NoseValues_Z, Method );
PhiVal = (-180:1:180);

NoseValues_XVal = NoseValues_XInterp(PhiVal);
NoseValues_YVal = NoseValues_YInterp(PhiVal);
NoseValues_ZVal = NoseValues_ZInterp(PhiVal);

       

% -------------- Plot --------------------
              
figure;
hold on
    % Plotting Values
%     scatter3( NoseValues_X, NoseValues_Y, NoseValues_Z, 50, PhiList, 'filled')
    NoseVal = scatter3( NoseValues_XVal, NoseValues_YVal, NoseValues_ZVal, 50, PhiVal, 'filled');
    axis equal
    
    % Plotting Seasons
    Offset = 0*2;
    Summer = scatter3( Offset+ NoseValues_XVal(PhiVal==0),  NoseValues_YVal(PhiVal==0),  NoseValues_ZVal(PhiVal==0), 300, 'filled');
    Winter = scatter3( Offset+ NoseValues_XVal(PhiVal==180),  NoseValues_YVal(PhiVal==180),  NoseValues_ZVal(PhiVal==180), 300, 'filled');
    Autumn = scatter3( Offset+ NoseValues_XVal(PhiVal==90),  NoseValues_YVal(PhiVal==90),  NoseValues_ZVal(PhiVal==90), 300, 'filled');
    Spring = scatter3( Offset+ NoseValues_XVal(PhiVal==-90),  NoseValues_YVal(PhiVal==-90),  NoseValues_ZVal(PhiVal==-90), 300, 'filled');
    
    ColourSummer = cbrewer2('Oranges', 20);
    set(Summer, 'MarkerEdgeColor', ColourSummer(15, :), 'LineWidth', 3, 'MarkerFaceColor', 'none');
    ColourSpring = cbrewer2('YlGn', 20);
    set([Spring, Autumn], 'MarkerEdgeColor', ColourSpring(15, :), 'LineWidth', 3, 'MarkerFaceColor', 'none');
    ColourWinter = cbrewer2('Greys', 20);
    set([Winter], 'MarkerEdgeColor', ColourWinter(15, :), 'LineWidth', 3, 'MarkerFaceColor', 'none');
    
    
    
    % Colourbar
    cb = colorbar;
    ColorMap = vertcat( cbrewer2( 'BuPu', 100 ), flipud(cbrewer2( 'BuPu', 100 )) );
    colormap(ColorMap)    
    set(cb, 'LineWidth', 1, 'FontName', 'Helvetica', 'XTick', [-180:45:180])

    % ---- Plot options
    axes = gca;
    set(axes, 'FontName'   , 'Helvetica',...
        'XAxisLocation', 'bottom',...
        'YAxisLocation', 'left',...
        'FontSize', 16,...
        'LineWidth'   , 1 ...
        );
    yticks([-6:2:6]);
    
    Ylabel = ylabel('Y (R_p)');
    Xlabel = xlabel('X (R_p)');
    Zlabel = zlabel('Z (R_p)');
    ColourBackground = [0, 0, 0, 0.3];
    set(gca,'Color',ColourBackground)
    set(gcf,'color','w');
    
    % ---- Axes and Grid
    axis([22, 29, -4, 4, -8, 8])
    grid on
       
    export_fig NoseLocus_XY.pdf -opengl -q101
    
    
%%  Flapping of Current Sheet    
        
	Mtilted = -horzcat( cos(PhiList*pi/180) .* sin(27*pi/180), sin(PhiList*pi/180) .* sin(27*pi/180), 0.*PhiList + cos(27*pi/180) );
	NoseValues = horzcat( NoseValues_X, NoseValues_Y, NoseValues_Z );
       
    LambdaSol = 0 .* PhiList;
    for kPhi = 1:length(PhiList)
        
        Mtilted_Phi = Mtilted(kPhi,:);
        NoseValues_Phi = NoseValues(kPhi,:);
    
        M_NormalLine = @(lambda) NoseValues_Phi + lambda .* Mtilted_Phi;
        funcOptim = @(lambda) sqrt( (M_NormalLine(lambda) * [1;0;0]).^2 + (M_NormalLine(lambda) * [0;1;0]).^2 + (M_NormalLine(lambda) * [0;0;1]).^2 );
    
        LambdaGuess = 10;    
        LambdaSol(kPhi) = fminsearch(funcOptim, LambdaGuess);
        
    end
    
    LambdaInterp = griddedInterpolant(PhiList, LambdaSol);
    LambdaVal = LambdaInterp(PhiVal) .' ;
    
    figure
    hold on
    scatter(PhiList, sqrt( ((LambdaSol.* Mtilted)* [1;0;0]).^2 + ((LambdaSol.* Mtilted)* [0;1;0]).^2 + ((LambdaSol.* Mtilted)* [0;0;1]).^2) );
    
    
% -- Position of Intersection w/t Sun-direction
    M_Sol = NoseValues + LambdaSol .* Mtilted;
    
    Method = 'spline';
    NoseInterp.X = griddedInterpolant(PhiList, NoseValues(:,1), Method );
    NoseInterp.Y = griddedInterpolant(PhiList, NoseValues(:,2), Method );
    NoseInterp.Z = griddedInterpolant(PhiList, NoseValues(:,3), Method );
    NoseInterp.XVal = NoseInterp.X(PhiVal);
    NoseInterp.YVal = NoseInterp.Y(PhiVal);
    NoseInterp.ZVal = NoseInterp.Z(PhiVal);
    
    MSolInterp.X = griddedInterpolant(PhiList, M_Sol(:,1), Method );
    MSolInterp.Y = griddedInterpolant(PhiList, M_Sol(:,2), Method );
    MSolInterp.Z = griddedInterpolant(PhiList, M_Sol(:,3), Method );
	MSolInterp.XVal = MSolInterp.X(PhiVal);
	MSolInterp.YVal = MSolInterp.Y(PhiVal);
	MSolInterp.ZVal = MSolInterp.Z(PhiVal);
    
    
    
    figure;
    hold on
        % Position Intersection
        Intersection = scatter3( MSolInterp.XVal, MSolInterp.YVal, MSolInterp.ZVal, 100, PhiVal );
        set(Intersection, 'Marker', '.');
        axis equal
        colorbar
        % Position Nose
        Nose = scatter3( NoseInterp.XVal, NoseInterp.YVal, NoseInterp.ZVal, 100, PhiVal, 'filled'  );
        set(Nose, 'Marker', 'o');

    % Colourbar
    cb = colorbar;
    ColorMap = vertcat( cbrewer2( 'BuPu', 100 ), flipud(cbrewer2( 'BuPu', 100 )) );
    colormap(ColorMap)    
    set(cb, 'LineWidth', 1, 'FontName', 'Helvetica', 'XTick', [-180:45:180])

    % ---- Plot options
    axes = gca;
    set(axes, 'FontName'   , 'Helvetica',...
        'XAxisLocation', 'bottom',...
        'YAxisLocation', 'left',...
        'FontSize', 16,...
        'LineWidth'   , 1 ...
        );
    yticks([-6:2:6]);
    
    Ylabel = ylabel('Y (R_p)');
    Xlabel = xlabel('X (R_p)');
    Zlabel = zlabel('Z (R_p)');
    ColourBackground = [0, 0, 0, 0.3];
    set(gca,'Color',ColourBackground)
    set(gcf,'color','w');
    
    % ---- Axes and Grid
    axis([21, 28, -4, 4, -14, 14])
    grid on
           
    export_fig NoseEquatorLocus.pdf -opengl -q101
    
    %{
    % TEST: Latitude of Sun, from Arridge2008
    Norm = sqrt( NoseValues_Dusk_Tot(:,1).^2 + NoseValues_Dusk_Tot(:,2).^2 + NoseValues_Dusk_Tot(:,3).^2 );
    tYr = linspace(0 ,1, length(Norm)) .';
    ThetaSun = -1.371 - 25.69 .* cos(2.816 + 0.213499 .* tYr) + 1.389 .* cos(5.4786 + 0.426998 .* tYr);
    X = 5.* cos(ThetaSun .* pi/180);
    Y = 5.* sin(ThetaSun .* pi/180);
    figure;
        scatter(X, Y)
        axis equal
    %}
    
    
    NM = NoseValues - M_Sol;
    NM_Sol_Norm = sqrt( NM(:,1) .^2 + NM(:,2) .^2 + NM(:,3) .^2 );
    figure;
    scatter3( NM(:,1), NM(:,2), NM(:,3), 100, PhiList, 'filled' )
    axis equal
    colorbar
%     axis([10 30 -10 10 -15 15])
    
    
    Cond1 = (-90<=PhiList & PhiList<=90);
    HN = M_Sol + NM_Sol_Norm .* Mtilted ./ sqrt( (Mtilted * [1;0;0]).^2 + (Mtilted * [0;1;0]).^2 + (Mtilted * [0;0;1]).^2 );
    for kPhi=1:length(PhiList)
        if (-90<=PhiList(kPhi) && PhiList(kPhi)<=90)
            HN( kPhi,3 ) = - HN( kPhi,3 );
        end
    end
    figure;
    scatter3( HN(:,1), HN(:,2), -HN(:,3), 100, PhiList, 'filled' )
    axis equal
    colorbar
    
              
%{           
              

% NoseValues_Dawn = [
%              [ 24.9037         0   -5.4299 ];...
%              [ 24.7479   -0.8942   -5.6892 ];...
%              [ 24.7092   -1.5327   -5.0235 ];...
%              [ 24.8583   -1.6180   -3.7839 ];...
%              [ 25.0487   -1.0586   -2.1180 ];...
%              [ 25.1740   -0.6790   -1.3326 ];...
%              [ 25.2882   -0.1974   -0.3948 ];...
%              [ 25.4528    0.6647    1.5544 ];...
%              [ 25.3974    1.0152    3.3276 ];...
%              [ 25.1932    0.7604    4.8377 ];...
%              [ 24.9037         0    5.4299 ];...
%              ];

% -------------- Dusk Side --------------------

PhiList_Dusk = horzcat( -20:-20: -80, -90).' ;




% NoseValues_Dusk =  [
%                     [25.1999, 0.7606,   -4.8390];...
%                     [25.3646    1.0139   -3.3233];...
%                     [25.4263    0.6640   -1.5528];...
%                     [25.2664   -0.1972    0.3945];...
%                     [25.1546   -0.6785    1.3316];...
%                     [25.0170   -1.0376    2.0759];...
%                     [24.8463   -1.6172    3.7820];...
%                     [24.7390   -1.5477    5.0727];...
%                     [24.7755   -0.8952    5.6955];...
%                     ];




% -------------- Concatenation -----------------
PhiList = vertcat(flipud(PhiList_Dusk), PhiList_Dawn);
NoseValues = vertcat(flipud(NoseValues_Dawn), NoseValues_Dawn);
NoseValues_X =  NoseValues(:, 1);
NoseValues_Y =  NoseValues(:, 2);
NoseValues_Z =  NoseValues(:, 3);

% -------------- Interpolation -----------------

NoseValues_XInterp = griddedInterpolant( PhiList, NoseValues_X, 'linear') ;
NoseValues_YInterp = griddedInterpolant( PhiList, NoseValues_Y, 'linear' ) ;
NoseValues_ZInterp = griddedInterpolant( PhiList, NoseValues_Z, 'linear' ) ;

PhiList_Interp = horzcat(-180:10:180) .' ;
NoseValues_XVal = NoseValues_XInterp(PhiList_Interp);
NoseValues_YVal = NoseValues_YInterp(PhiList_Interp);
NoseValues_ZVal = NoseValues_ZInterp(PhiList_Interp);


Cond = ( -180 <= PhiList & PhiList <= 180) ;

figure;
hold on
    scatter3( NoseValues_X(Cond), NoseValues_Y(Cond), NoseValues_Z(Cond), 50, PhiList(Cond), 'filled')
    scatter3( NoseValues_XVal,  NoseValues_YVal,  NoseValues_ZVal, 50, PhiList_Interp, 'filled')
    Summer = scatter3( NoseValues_XVal(PhiList_Interp==0),  NoseValues_YVal(PhiList_Interp==0),  NoseValues_ZVal(PhiList_Interp==0), 200, 'filled');
    Winter = scatter3( NoseValues_XVal(PhiList_Interp==180),  NoseValues_YVal(PhiList_Interp==180),  NoseValues_ZVal(PhiList_Interp==180), 200, 'filled');
    Autumn = scatter3( NoseValues_XVal(PhiList_Interp==90),  NoseValues_YVal(PhiList_Interp==90),  NoseValues_ZVal(PhiList_Interp==90), 200, 'filled');
    Spring = scatter3( NoseValues_XVal(PhiList_Interp==-90),  NoseValues_YVal(PhiList_Interp==-90),  NoseValues_ZVal(PhiList_Interp==-90), 200, 'filled');
    cb = colorbar;
    C = 7;     
    axis equal
    ColorMap = vertcat( cbrewer2( 'BuPu', 100 ), flipud(cbrewer2( 'BuPu', 100 )) );
    colormap(ColorMap)    
    
        % ---- Plot options
        axes = gca;
        set(axes, 'FontName'   , 'Helvetica',...
            'XAxisLocation', 'bottom',...
            'YAxisLocation', 'left',...
            'FontSize', 16,...
            'LineWidth'   , 1 ...
            );
        Ylabel = ylabel('Y (R_p)');
        Xlabel = xlabel('X (R_p)');
        Zlabel = zlabel('Z (R_p)');
        ColourBackground = [0, 0, 0, 0.3];
        set(gca,'Color',ColourBackground)
        set(gcf,'color','w');
        
%         set(cb, 'LineWidth', 1, 'FontName', 'Helvetica', 'XTick', [-180:45:180])

    axis([22, 27, -5, 5, -7, 7])
    grid on




%{


PhiList = horzcat( 0:20:80, 90, 100:20:260, 270, 280:20:360  ).' ;
Nose_List = 0.* PhiList;


NoseValues_Psw0_01 = [
             [ 24.9037         0   -5.4299 ];...
             [ 24.7479   -0.8942   -5.6892 ];...
             [ 24.7092   -1.5327   -5.0235 ];...
             [ 24.8583   -1.6180   -3.7839 ];...
             [ 25.0487   -1.0586   -2.1180 ];...
             [ 25.1740   -0.6790   -1.3326 ];...
             [ 25.2882   -0.1974   -0.3948 ];...
             [ 25.4528    0.6647    1.5544 ];...
             [ 25.3974    1.0152    3.3276 ];...
             [ 25.1932    0.7604    4.8377 ];...
             [ 24.9037         0    5.4299 ];...
             ];
         
PhiList_Dawn = horzcat( -20:-20:-60 ).' ;
NoseValues_Psw0_01_Dawn = [
                          [25.1999, 0.7606,   -4.8390];...
                          [25.3646    1.0139   -3.3233];...
                          [25.4263    0.6640   -1.5528];...
                          ];

NoseValues_rMP25 = [ [ 24.4212;  0;  -5.3694], ...
               [ 24.3658; -0.8735; -5.5571 ], ...
               [ 24.4695; -1.4918; -4.8894 ], ...
               [ 24.6573; -1.6049; -3.7532 ], ...
               [ 24.8633; -1.0312; -2.0631 ], ...
               [ 24.9559; -0.6533; -1.2821], ...
               [ 24.9941; -0.1755; -0.3512 ],...
               [ 24.9340; 0.7027; 1.6434 ],...
               [ 24.7342; 1.0530;  3.4513 ],...
               [24.5343; 0.7405; 4.7111 ],...
               [24.4212; 0; 5.3694 ],...
               ].';
           
NoseValues = vertcat(NoseValues_Psw0_01_Dawn, NoseValues_Psw0_01)   ;
PhiList = vertcat(PhiList_Dawn, horzcat( 0:20:80, 90, 100:20:180).' );
           
           
%  % Forcing DD Symmetry
%  NoseValues = NoseValues_Psw0_01 ;  
%  NoseValues_Dusk = NoseValues;  
%  NoseValues_Dusk(:, 2) = -NoseValues_Dusk(:, 2);
%  NoseValues = vertcat(NoseValues, flipud(NoseValues_Dusk(1:end-1, :)));          
%          
%                [ 24.4588; -0.8768; 5.5783], ...
%                [ 24.6841; -1.4917; 4.8893], ...
%                [ 25.0476; -1.5950; 3.7301], ...
%                [ 25.4100; -1.0339; 2.0685], ...
%                [ 25.6010; -0.6498; 1.2753], ...
%                [ 25.7605; -0.1809; 0.3620], ...
%                [ 25.9370; 0.6058; -1.4168], ...
%                [ 25.8913; 0.9411; -3.0846], ...
%                [ 25.6329; 0.7018; -4.4649], ...
%                [ 24.4212;  0;  -5.3694], ...
%                ];

Noseinterp.XInterp = griddedInterpolant(PhiList, NoseValues(:, 1), 'spline');
Noseinterp.YInterp = griddedInterpolant(PhiList, NoseValues(:, 2), 'spline');
Noseinterp.ZInterp = griddedInterpolant(PhiList, NoseValues(:, 3), 'spline');

PhiList_2 = (0:1:360);
Noseinterp.X = Noseinterp.XInterp(PhiList_2);
Noseinterp.Y = Noseinterp.YInterp(PhiList_2);
Noseinterp.Z = Noseinterp.ZInterp(PhiList_2);

figure;
hold on
    scatter3(Noseinterp.X, Noseinterp.Y, Noseinterp.Z , 200, PhiList_2, 'filled')
    scatter3( NoseValues(:, 1),NoseValues(:, 2),  NoseValues(:, 3), 200, PhiList, 'filled')
    cb = colorbar;
    C = 7;     
    axis equal
    ColorMap = vertcat( cbrewer2( 'BuPu', 100 ), flipud(cbrewer2( 'BuPu', 100 )) );
    colormap(ColorMap)    
    
        % ---- Plot options
        axes = gca;
        set(axes, 'FontName'   , 'Helvetica',...
            'XAxisLocation', 'bottom',...
            'YAxisLocation', 'left',...
            'FontSize', 16,...
            'LineWidth'   , 1 ...
            );
        Ylabel = ylabel('Y (R_p)');
        Xlabel = xlabel('X (R_p)');
        Zlabel = zlabel('Z (R_p)');
        ColourBackground = [0, 0, 0, 0.3];
        set(gca,'Color',ColourBackground)
        set(gcf,'color','w');
        
        set(cb, 'LineWidth', 1, 'FontName', 'Helvetica', 'XTick', [0:45:360])

    axis([22, 27, -5, 5, -7, 7])
    
    export_fig NoseLocus_XOZ.pdf -opengl -q101
    
    
    
    
    
    
    
    
figure;
hold on
    for kPhi = 1:length(PhiList)        
        Nose = NoseValues(:, kPhi);
        scatter3(Nose(1), Nose(2), Nose(3), 100, 'filled')
    end
axis equal