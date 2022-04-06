
function PhiRmpPsw_Interp = PhiRmpPsw_Interp( Conditions_PhiTilt_Psw_rMP )

    PhiList = (0:20:180).' ;
    rMPList = (15:2:35) .' ;

    PswValues = [...
                [ 0.1263; 0.0619; 0.03534; 0.0226; 0.0155; 0.0108; 0.0080; 0.0061; 0.004831; 0.0039; 0.003235 ], ...
                [ 0.1253; 0.0612; 0.0342; 0.0217;  0.0151; 0.0106; 0.0078; 0.0059; 0.0046; 0.003775; 0.003130 ],...
                [ 0.1217; 0.0620; 0.0333; 0.0207;  0.0145; 0.0104; 0.0075; 0.0057; 0.00446; 0.00362; 0.003015 ],...
                [ 0.1239; 0.0631; 0.0338; 0.0201; 0.0141; 0.0104; 0.0074; 0.0056; 0.00435; 0.00353; 0.00296], ...
                [ 0.1189; 0.0652; 0.0356; 0.0209; 0.0142; 0.01034; 0.00762; 0.0058; 0.00444; 0.00358; 0.00298], ...
                [ 0.1169; 0.0669; 0.0368; 0.0217; 0.0146; 0.0105; 0.00783; 0.0059; 0.0046; 0.0037; 0.00306] ,...
                [ 0.1240; 0.0675; 0.0372; 0.0220; 0.0152; 0.0109; 0.0080; 0.00606; 0.00474;  0.00383; 0.00316], ...
                [ 0.1304; 0.0660; 0.0368; 0.0228; 0.0157; 0.01115; 0.00814; 0.0062; 0.00488; 0.00392; 0.00324 ], ...
                [ 0.1306; 0.0639; 0.0366; 0.0230; 0.0158; 0.0111; 0.00813; 0.00623; 0.0049; 0.00396; 0.00327 ], ...
                [ 0.1263; 0.0619; 0.0354; 0.0226; 0.0155; 0.0108; 0.008; 0.0061; 0.00482; 0.0039; 0.003235 ], ...
                ];

    [PhiGrid, rMPGrid] = ndgrid(PhiList, rMPList);
    log10Psw_Interp = griddedInterpolant( PhiGrid, log10(rMPGrid), log10(PswValues.'), 'spline');
    
    if isnan(Conditions_PhiTilt_Psw_rMP(2))
        
        Phi = Conditions_PhiTilt_Psw_rMP(1);
        logRmp = log10(Conditions_PhiTilt_Psw_rMP(3));
        log10Psw = log10Psw_Interp(Phi, logRmp);
        
        PhiRmpPsw_Interp = 10^log10Psw;
        
    elseif isnan(Conditions_PhiTilt_Psw_rMP(3))
        
        Phi = Conditions_PhiTilt_Psw_rMP(1);
        logPsw = log10( Conditions_PhiTilt_Psw_rMP(2) );
        
        Cont = contourc( PhiList, log10(rMPList), log10(PswValues), [logPsw logPsw]);
        logrMPInterp = griddedInterpolant( Cont(1,2:end), Cont(2,2:end) );
        logRmp = logrMPInterp(Phi);
        
        PhiRmpPsw_Interp = 10^logRmp;
        
    end
    

end


%{
PhiTilt = 40;
rMP = 35;
PswGuess = 0.00296;

FuncRmp = @(Psw) testPsw(PhiTilt, rMP, Psw);
PswSol = fsolve( FuncRmp, PswGuess )








figure;
hold on
for kPhi = 1:length(PhiList)
    Phi = PhiList(kPhi);
    PswGrid()
    
	X = log10( PswValues(:, kPhi) );
	Y = log10(15:2:35).' ;
	scatter( X, Y, 100, 'filled')
end
        XInterp =  horzcat(1+0.*X, X);
        Beta = XInterp \ Y;
        YInterp = XInterp * Beta;
        plot( X, YInterp)
%}
        
        
        