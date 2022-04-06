

%% ON THE SIDE

%{

%% ------------------------------------------------------------------------
% Finding Position of Anchor Points: Stagnation Point and Terminator Point

% ----  Finding B.v=0 locus in NMM meridian plane + Stagnation Point
%     ParaSystem = StagnationPoint(ParaSystem);
    ParaSystem = NoseDetermination(ParaSystem);
   
% ----  Finding Bxv=0 locus in NMM meridian plane + Terminator Point
    ParaSystem = TerminatorPoint(ParaSystem);



%% -------------------------------------------------------------------------------
% Solving PB in NMM: Piece-Wise Approach || Stagnation Point Centered!

% ----  Branch 1: From Stagnation Point, Downstream


% From SSN: first branch

CoordinateDetails.Type = 'Spherical_MB';
CoordinateDetails.Units = 'r0';

% Finding Value of drdtheta at Stagnation Point 
ex = [1; 0; 0];
StPoint_r0 = ParaSystem.Nose.KSMU_Rp * ParaSystem.Rp / ParaSystem.r0;
StPointDirection = StPoint_r0 ./ norm(StPoint_r0);
AlphaStPoint = acos(dot(ex, StPointDirection));
drdtheta_StPoint = -norm(StPoint_r0) * tan(AlphaStPoint);


% Defining PB in NMM Plane 

drdphi_NMM = 0;
drdtheta_ini = drdtheta_StPoint;
PB_NMM = @(theta, r, dr) PressureBalance_Equation_StPoint(CoordinateDetails, theta, r, dr, drdphi_NMM, ParaSystem);
PB_NMM(0, norm(StPoint_r0), drdtheta_ini)

drFunc = @(dr) PB_NMM(0, norm(StPoint_r0), dr);
drFunc(drdtheta_ini)
fsolve(drFunc, 0)

% Determining First Branch  

ThetaSpan = (0:90).'.*DegToRad;
r_ini = norm(StPoint_r0);
dr_ini = drdtheta_ini;

[theta, r] = ode15i(PB_NMM, ThetaSpan, r_ini, dr_ini);

rVector = r * ParaSystem.r0 / ParaSystem.Rp;
PhiVector = 0*theta + pi/2-ParaSystem.AlphaDipole;
ThetaVector = theta-AlphaStPoint;
[Y, Z, X] = sph2cart(PhiVector, pi/2-ThetaVector, rVector);

% Test : Explicit
func = @(theta, r) drMeridianCAN_1(CoordinateDetails, theta, r, ParaSystem);
rGuess = r_ini;
[Theta, rSol] = ode45(func, ThetaSpan, rGuess)


%% TEST

hold on
NMMBranch1 = plot3(X, Y, Z);
set(NMMBranch1, 'linewidth', 2)

% ----  Branch 2: From Terminator, Upstream
ThetaSpan = (90:-30).'.*DegToRad;
StPoint_A1 = ParaSystem.A1.TiltedDipole.KSM_Rp * ParaSystem.Rp / ParaSystem.r0;

r_ini = norm(StPoint_A1);
drFunc = @(dr) PB_NMM(pi/2, r_ini, dr);
drGuess = 10;

figure;
hold on
drList = -10:10;
for k =1:length(drList)
    scatter(drList(k), drFunc(drList(k)))

end

drFunc(drGuess)
drSol = fsolve( drFunc, drGuess)

dr_ini = drdtheta_ini;

[theta, r] = ode15i(PB_NMM, ThetaSpan, r_ini, dr_ini);

% figure;
% plot(ThetaVector * 180/pi, X)
% 
% ThetaVector(X==max(X)) * 180/pi
% AlphaStPoint * 180/pi


%% -------------------------------------------------------------------------------
% Spectral Method

dPhi = 4*DegToRad;
L = pi;
Phi = (-L/2:dPhi:L/2-dPhi).';
r0 = 0*Phi + norm(ParaSystem.Nose.KSMU_Rp * ParaSystem.Rp / ParaSystem.r0);

ThetaCusp = 71;

dTheta = 1*DegToRad;
ThetaMax = 10;
NTheta = ThetaMax / (dTheta*RadToDeg);

Thetak = 0;
Data = zeros( length(Phi),  floor(NTheta));
Data(:, 1) = r0;

for k = 1:NTheta
    
    display(num2str(floor((k/NTheta) * 100)))
%     Thetak = (k)*dTheta;
    PB_Func = @(Theta, r) PB_Spectral( Theta, r, L, Phi, ThetaCusp, ParaSystem );
    [Theta, rSol] = ode45( PB_Func, [Thetak, Thetak+dTheta], r0 );
    Thetak = Thetak + dTheta;
    
    r0 = real(rSol(end, :));
    Data(:, k+1) = r0;
end


    Phi = vertcat(Phi, Phi(1));
    Data = vertcat(Data, Data(1,:));
    
    figure;
    ThetaSpan = (0:dTheta:ThetaMax* DegToRad).'
    [PhiGrid, ThetaGrid] = ndgrid(Phi, ThetaSpan);
    [Y, Z, X] = sph2cart(PhiGrid, pi/2-ThetaGrid, Data);
    Surface1 = surf(-X, Y, Z, -X);
    axis equal
%     set(Surface1, 'FaceColor', (Color_NMM(30,:)), 'FaceAlpha', 0.5)






%% TEST: Compression guess surface
% C = 0.7;
% Compression = (C-1)*abs(sin(ParaGrid.PhiGrid)) + 1;
% rGrid = Compression.*rGrid;

%% TEST

%{
[X, Y, Z] = sph2cart(ParaGrid.PhiGrid, pi/2-ParaGrid.ThetaGrid, rGrid);
figure;
surf(X, Y, Z)
axis equal


ResidualGrid = PB_Grid(rGrid, ParaGrid, ParaSystem);
figure;
pcolor(ResidualGrid)
colorbar
%}

% TEST
% rGrid = rGrid(:, 2:end);
% 
% ParaGrid.PhiList = (0: ParaGrid.DeltaPhi * RadToDeg :360).' .* DegToRad;
% ParaGrid.ThetaList = (ParaGrid.DeltaTheta * RadToDeg: ParaGrid.DeltaTheta * RadToDeg :90).' .* DegToRad;
% 
% [ParaGrid.PhiGrid, ParaGrid.ThetaGrid] = ndgrid(ParaGrid.PhiList, ParaGrid.ThetaList);
% TEST









MaxIter = 1500;
MaxEval = 5e5;
ParaGrid.epsilon = 1e-5;
PBGrid_Func = @(r) PB_Grid(r, ParaGrid, ParaSystem);
Algorithms = ["levenberg-marquardt"; "Trust-Region-Dogleg"];
opts = optimoptions('fsolve', 'Display','iter', 'Maxiter', MaxIter, 'MaxFunctionEvaluations', MaxEval, 'SpecifyObjectiveGradient', true, 'CheckGradients', false, 'Algorithm', Algorithms(2));
% [f, J] = PBGrid_Func(rGrid)
% 
% figure;
% error = abs((J-grad_fd));
% pcolor(error); shading flat
% colorbar

rSol = fsolve(PBGrid_Func, rGrid, opts)
rGrid = rSol;
% Test
% rSol(:,1) = norm(Nose_r0);
% PhiList = ParaGrid.PhiList;
% Position180 = find(PhiList == pi);
% Position90 = (Position180+1)/2;
% Position270 = find(PhiList == 3*pi/2);
% 
% 
% rSol(Position90+1:Position180,:) = flipud(rSol(1:Position90-1,:));
% rSol(Position270+1:end,:) =  flipud(rSol(Position180:Position270-1,:));
% % rSol(Position90, :) = (1/2) .* (rSol(Position90-1, :) + rSol(Position90+1, :) );
% % rSol(Position270, :) = (1/2) .* (rSol(Position270-1, :) + rSol(Position270+1, :) );
% Test







% figure;
% surf(ZNose, XNose, YNose)
% axis equal

%{
YNoseInter = griddedInterpolant(ParaGrid.PhiGrid, ParaGrid.ThetaGrid, YNose);
ZNoseInter = griddedInterpolant(ParaGrid.PhiGrid, ParaGrid.ThetaGrid, ZNose);
XNoseInter = griddedInterpolant(ParaGrid.PhiGrid, ParaGrid.ThetaGrid, XNose);

ThetaVec = (0:1:90) .* DegToRad;
PhiVec = (0:1:360) .* DegToRad;
[PhiGrid, ThetaGrid] = ndgrid(PhiVec, ThetaVec);
YNose = YNoseInter(PhiGrid, ThetaGrid);
ZNose = ZNoseInter(PhiGrid, ThetaGrid);
XNose = XNoseInter(PhiGrid, ThetaGrid);
%}








%% In PB_Grid.m

% TEST
% % 
% PhiList = ParaGrid.PhiList;
% Position180 = find(PhiList == pi);
% Position90 = (1/2)*(Position180+1);
% Position270 = find(PhiList == 3*pi/2);
% 
% % 
% rGrid(Position90+1:Position180,:) = flipud(rGrid(1:Position90-1,:));
% rGrid(Position270+1:end,:) =  flipud(rGrid(Position180:Position270-1,:));
% rGrid(Position90, :) = (1/2) .* (rGrid(Position90-1, :) + rGrid(Position90+1, :) );
% rGrid(Position270, :) = (1/2) .* (rGrid(Position270-1, :) + rGrid(Position270+1, :) );

% % 


% TEST




% rGrid_pNphi = horzcat( rGrid_pNphi, 2*rGrid_pNphi(:, end) - rGrid_pNphi(:, end-1));
% rGrid_pNphi(:, end) = 2*rGrid(:, end) - rGrid_mNphi(:, end);
% drdthetaGrid = (1/2) .* (drdthetaGrid_BWD + drdthetaGrid_FWD);
% drdthetaGrid = drdthetaGrid_Central;
% drdthetaGrid(:,1) = drdthetaGrid_FWD(:,1);
% drdthetaGrid(:,2:end) = drdthetaGrid_Central(:,2:end);
drdthetaGrid = drdthetaGrid_FWD;



% TEST
% rGrid_p1 = rGrid(2:end, :);
% rGrid_p1 = vertcat(rGrid_p1, 2*rGrid(end, :) - rGrid(end-1, :) );
% 
% rGrid_m1 = rGrid(1:end-1, :);
% rGrid_m1 = vertcat( 2*rGrid(1, :) - rGrid(2, :), rGrid_m1 );
% TEST




% TEST: Introduce Phi symmetry
% PhiList = ParaGrid.PhiList;
% Position90 = find(PhiList == pi/2);
% Position180 = find(PhiList == pi);
% Position270 = find(PhiList == 3*pi/2);
% 
% nGrid_ePhi(Position90,:) = 0;
% nGrid_ePhi(Position270,:) = 0;

% nGrid_ePhi(Position90+1:Position180,:) = nGrid_ePhi(1:Position90-1,:);
% nGrid_ePhi(Position270+1:end,:) = nGrid_ePhi(Position180:Position270-1,:);
% 
% TEST: Introduce Phi symmetry






% Conversion of nGrid from Nose_Spherical to Nose_Cartesian
% nX_Nose = nGrid_er .* sin(ThetaGrid) .* cos(PhiGrid) +  nGrid_eTheta .* cos(ThetaGrid) .* cos(PhiGrid) +  nGrid_ePhi .* (-sin(PhiGrid));
% nY_Nose = nGrid_er .* sin(ThetaGrid) .* sin(PhiGrid) +  nGrid_eTheta .* cos(ThetaGrid) .* sin(PhiGrid) +  nGrid_ePhi .* (cos(PhiGrid));
% nZ_Nose = nGrid_er .* cos(ThetaGrid) +  nGrid_eTheta .* (-sin(ThetaGrid));

% nX_Nose = nZ_Nose;
% nY_Nose = nX_Nose;
% nZ_Nose = nY_Nose;










%% ------------------------------------------------------------------------
% Script

 
%{    
%% -------------------------------------------------------------------------------
% Analytical Position of Nose

    PhiList = (0:5:360).' .* DegToRad;
    rList = 0*PhiList;
    ThetaList =  0*PhiList;
    ResidualList = 0*PhiList;

    for kPhi = 1:length(PhiList)

        Phi = PhiList(kPhi);

        ParaSystem.Tilt.Phi = Phi;
        ParaSystem.Tilt.RotationPhi = [ cos(ParaSystem.Tilt.Phi) -sin(ParaSystem.Tilt.Phi) 0 ;...
                                    sin(ParaSystem.Tilt.Phi)  cos(ParaSystem.Tilt.Phi) 0 ;...
                                    0              0             1];

        Theta = ParaSystem.Tilt.Theta;
        TanTheta = sec(Phi) * ( cot(Theta) - sqrt( (cos(Phi))^2+(cot(Theta))^2 )  )  ;
        NoseTheta = atan(TanTheta);
        ThetaList(kPhi) = NoseTheta;

        fFunction = cos(Theta)*(3*(sin(NoseTheta))^2-1)+(3)*cos(Phi)*sin(Theta)*sin(NoseTheta)*cos(NoseTheta);
        rSol = (2*ParaSystem.M)^(1/3) * ( (sin(Phi))^2*(sin(Theta))^2 + fFunction^2 ) ^ (1/6);

        PBNose = @(r) PB_Nose(r, NoseTheta, ParaSystem);
        rList(kPhi) = rSol;
        ResidualList(kPhi) = PBNose(rSol);
    end

    PhiNMM = acos(sin(PhiList) .* sin(Theta));

    rList_Rp = rList .* ParaSystem.r0 ./ ParaSystem.Rp;
    XList = rList_Rp .* cos(ThetaList);
    YList = rList_Rp .* sin(ThetaList) .* cos(PhiNMM);
    ZList = rList_Rp .* sin(ThetaList) .* sin(PhiNMM);

    %{
    hold on
        NoseLocus = scatter3(XList, YList, ZList, 40, PhiList * RadToDeg, 'filled');
        Color_NMMPlane = cbrewer2('Greens', 40);
        Opacity = 0.5;
        set(NoseLocus, 'MarkerFaceColor', Color_NMMPlane(10,:), 'MarkerFaceAlpha', Opacity)
    %}

% Resetting Rotation Angles
    ParaSystemRef = System_ini(Planet, Inclination);
    
    ParaSystem.Tilt.Theta =  ParaSystemRef.Tilt.Theta; % Dipole tilt, in rad
    ParaSystem.Tilt.Phi =  ParaSystemRef.Tilt.Phi; % Dipole tilt, in rad
    ParaSystem.Tilt.RotationPhi = ParaSystemRef.Tilt.RotationPhi;
%}    
   





%}


%{ 

    
    
%{

% SCRIPT
%% -------------------------------------------------------------------------------
% Optimising Towards Pressure Balance: Setting Up Grid


    ParaGrid.DeltaPhi = 2.* DegToRad;
    ParaGrid.DeltaTheta = 4.* DegToRad;

%     ParaGrid.PhiList = (0: ParaGrid.DeltaPhi * RadToDeg :360).' .* DegToRad;
    ParaGrid.ThetaList = (0: ParaGrid.DeltaTheta * RadToDeg :90).' .* DegToRad;
    
    ParaGrid.PhiList = (-90: ParaGrid.DeltaPhi * RadToDeg : 90).' .* DegToRad;

    
    [ParaGrid.PhiGrid, ParaGrid.ThetaGrid] = ndgrid(ParaGrid.PhiList, ParaGrid.ThetaList);

        rGrid_ini = norm(Nose_r0) * (2./(1+cos(ParaGrid.ThetaGrid))).^(0.6);
%}
    
    
    
    
%{  
    %% TEST: Explicit Intergration in NMM
    
    
PBNMM1 = @(theta, r) drdtheta1(theta, r, ParaSystem);
PBNMM2 = @(theta, r) drdtheta2(theta, r, ParaSystem);

% First Branch: From SNN
r0 = norm(ParaSystem.Nose.KSMU_Rp .* ParaSystem.Rp ./ ParaSystem.r0);
ThetaSpan = (0.1:0.1:100).* DegToRad;

options = odeset('RelTol',1e-5, 'AbsTol',1e-5);
[ThetaSol, rSol] = ode45(PBNMM1, ThetaSpan, r0, options);
[X, Y] = pol2cart(ThetaSol, rSol);

figure;
plot(Y, X)
axis equal
axis([ 0 1.3 0 1.3])    

max( ThetaSol( rSol == real(rSol) ) ) .* 180/pi;
    

% Second Branch: From A1
ParaSystemNorth = System_ini(Planet, Inclination);
ParaSystemNorth = TerminatorPoint(ParaSystemNorth);
 
rA1 = norm(ParaSystemNorth.A1.TiltedDipole.KSM_Rp .* ParaSystem.Rp ./ ParaSystem.r0);
Thetar0 = acos( ParaSystem.Nose.KSMU_Rp(1) .* ParaSystem.Rp ./ ParaSystem.r0 ./ r0) .* RadToDeg;

ThetaA1 = acos( dot(ParaSystem.NoseDirection, ParaSystem.A1.TiltedDipole.KSM_Rp.* ParaSystem.Rp ./ ParaSystem.r0) ./ rA1)   .* RadToDeg;
ThetaSpan2 = (ThetaA1 : -0.01 : 50).' .* DegToRad;

options = odeset('RelTol',1e-7, 'AbsTol',1e-7);
% ThetaSol = fsolve( @(theta) PBNMM2( theta, rA1 ), ThetaA1 .* DegToRad, options)
% ThetaSpan2 = (ThetaSol .* RadToDeg : -0.01 : 50).' .* DegToRad;

[ThetaSol2, rSol2] = ode45(PBNMM2, ThetaSpan2, rA1, options);
% 
% PBNMM1( ThetaSpan2(1), rA1 )
% PBNMM2( ThetaSpan2(1), rA1 )
% 
% 
% ThetaSol = fsolve( @(theta) PBNMM2( theta, rA1 ), ThetaA1 .* DegToRad, options)
% ThetaSpan2 = (ThetaSol .* RadToDeg : -0.01 : 50).' .* DegToRad;

% PBNMM2( ThetaSpan2(1), rA1 )
% PBNMM2( ThetaA1 .* DegToRad, rA1 )
% PBNMM2( (ThetaA1 + 0.01).* DegToRad, rA1 )
% 
% figure;
% thetaSpan = (90-3*0.01:0.01:90+3*0.01).* DegToRad;
% plot( thetaSpan .* RadToDeg, PBNMM2( thetaSpan, rA1 ) )


[X2, Y2] = pol2cart((ThetaSol2), (rSol2));

% figure;
hold on
plot(Y2, X2)
axis equal



% Test: finding shooting angle from A2

    ParaGrid.PhiList = pi/2;
    ParaGrid.ThetaList = ThetaA1 .* DegToRad;
    [ParaGrid.PhiGrid, ParaGrid.ThetaGrid] = ndgrid(ParaGrid.PhiList, ParaGrid.ThetaList);
    
    PBGrid_Func = @(r, drdtheta) PB_Grid_A2(r, drdtheta, ParaGrid, ParaSystem);
    fsolve( @(drdtheta) PBGrid_Func( rA1, drdtheta), 0);
    
    drdthetaSpan = (-2:0.01:2);
    rSpan = (rA1*0.8 :rA1*0.4/100 : rA1*1.2);
    [drGrid, rGrid] = ndgrid(drdthetaSpan, rSpan);
    PBGrid = PBGrid_Func( rGrid, drGrid);
    
    figure;
    hold on
    surf(rGrid, drGrid, PBGrid ./ max(PBGrid(:)))
    shading flat
    
    
    
    
    
    for kdr = 1:length(drdthetaSpan)
        for kr = 1:length(rSpan)
            drdtheta = drdthetaSpan(kdr);
            r = rSpan(kr);
            scatter3( drdtheta, r,  PBGrid_Func( r, drdtheta))
        end
    end
    
    %% TEST: Explicit Intergration in NMM

%}    


%}

%{ 
% PB in NMM, Northern Branch

%     ParaGrid.ThetaList = (0 + ParaGrid.DeltaTheta* RadToDeg  : ParaGrid.DeltaTheta * RadToDeg : ( (Nose_Theta+A1Theta) - ParaGrid.DeltaTheta) * RadToDeg).' .* DegToRad;
%     ParaGrid.ThetaList = (0 + ParaGrid.DeltaTheta* RadToDeg  : ParaGrid.DeltaTheta * RadToDeg : 70).' .* DegToRad;
%     ParaGrid.ThetaList = ((A1Theta-Nose_Theta) * RadToDeg : -ParaGrid.DeltaTheta * RadToDeg : 0).' .* DegToRad;


%}


%{

% PB in NMM, Northern Branch

    
    I1 = (1:length(ParaGrid.PhiList)*length(ParaGrid.ThetaList)).';
    J1 = I1;
    V1 = 0*I1 + 1;
    
    I2 = I1(:, 2:end);
    J2 = I2-1;
    V2 = 0*I2 + 1;
    
    I3 = I1(:, 1:end-1);
    J3 = I3+1;
    V3 = 0*I3+1;
        
    I = vertcat(I1, I2, I3);
    J = vertcat(J1, J2, J3);
    V = vertcat(V1, V2, V3);

    Jstr = sparse(I, J, V);


%}


%{ 

% PB in NMM, Northern Branch
   %%
   
   
%     rSolNMM_1 = vertcat(norm(Nose_r0), rSolNMM_1);
%     rNMM_Interp = griddedInterpolant( vertcat(0, ParaGrid.ThetaList), rSolNMM_1);
    

%     rNMM_Interp = griddedInterpolant(ParaGrid.ThetaList, rGrid);
    
%{
    ParaGrid.DeltaPhi = 2.* DegToRad;
    ParaGrid.DeltaTheta = 4.* DegToRad;

%     ParaGrid.PhiList = (0: ParaGrid.DeltaPhi * RadToDeg :360).' .* DegToRad;
%     ParaGrid.ThetaList = (0: ParaGrid.DeltaTheta * RadToDeg :90).' .* DegToRad;
    ParaGrid.ThetaList = (ParaGrid.DeltaTheta * RadToDeg: ParaGrid.DeltaTheta * RadToDeg :95).' .* DegToRad;


    ParaGrid.PhiList = (-90: ParaGrid.DeltaPhi * RadToDeg : 90).' .* DegToRad;

    
    [ParaGrid.PhiGrid, ParaGrid.ThetaGrid] = ndgrid(ParaGrid.PhiList, ParaGrid.ThetaList);
    
    rGrid = rNMM_Interp(ParaGrid.ThetaGrid);
  %}
        %{
    
%% TEST: NMM, branch 2
    
    ParaGrid.PhiList = -pi/2;
    ParaGrid.ThetaList = (0: ParaGrid.DeltaTheta * RadToDeg :90).' .* DegToRad;
    
    [ParaGrid.PhiGrid, ParaGrid.ThetaGrid] = ndgrid(ParaGrid.PhiList, ParaGrid.ThetaList);
    
    rGrid = norm(Nose_r0)* (2./(1+cos(ParaGrid.ThetaGrid))).^(kappaSol).';
    MaxIter = 30;
    MaxEval = 5e5;
    ParaGrid.epsilon = 1e-10;
    PBGrid_Func = @(r) PB_Grid_NMM(r, ParaGrid, ParaSystem);
    Algorithms = ["levenberg-marquardt"; "Trust-Region-Dogleg"; "Trust-Region" ];
    opts = optimoptions('fsolve', 'Display','iter', 'Maxiter', MaxIter, 'MaxFunctionEvaluations', MaxEval, 'SpecifyObjectiveGradient', false, 'CheckGradients', false, 'FiniteDifferenceType', 'central', 'Algorithm', Algorithms(1));
    
    rGuess = rGrid(:);
    rSolNMM_2 = fsolve(PBGrid_Func, rGuess, opts);
    rGridNMM_2 = rSolNMM_2;
    [X, Y] = pol2cart(ParaGrid.ThetaList, rSolNMM_2);
    figure;
    plot(Y, X)
    axis equal
        
%% Building Guess Surface

    rSolNMM_Interp_1 = griddedInterpolant(ParaGrid.ThetaList, rSolNMM_1);
    rSolNMM_Interp_2 = griddedInterpolant(ParaGrid.ThetaList, rSolNMM_2);
    
    ParaGrid.DeltaPhi = 5.* DegToRad;
    ParaGrid.DeltaTheta = 5.* DegToRad;
    ParaGrid.PhiList = (-90: ParaGrid.DeltaPhi * RadToDeg : 90).' .* DegToRad;
    ParaGrid.ThetaList = (0: ParaGrid.DeltaTheta * RadToDeg :90).' .* DegToRad;
    [ParaGrid.PhiGrid, ParaGrid.ThetaGrid] = ndgrid(ParaGrid.PhiList, ParaGrid.ThetaList);
        
    Flattening = 0.7;
    b = abs(rSolNMM_Interp_1(ParaGrid.ThetaList).*sin(ParaGrid.ThetaList) + rSolNMM_Interp_2(ParaGrid.ThetaList) .* sin(ParaGrid.ThetaList)) ./2;
    a = b ./ Flattening;
    e = sqrt(1-(b./a).^2);
    e(1) = e(2);
    

    rGrid = 0 .* ParaGrid.PhiGrid;
    
%     for kPhi = 1:length(ParaGrid.PhiList)
%         Phik = ParaGrid.PhiList(kPhi);
%         rGrid(kPhi, :) = (  (b) ./ sqrt(1-(e.*cos(Phik)).^2) ) ./ sin(ParaGrid.ThetaList) ;
%     end
%     rGrid (:, 1) = norm(ParaSystem.Nose.KSMU_Rp .* (ParaSystem.Rp ./ (ParaSystem.r0)));
%     rSol = rGrid;
    
    for kTheta = 1:length(ParaGrid.ThetaList)
            Thetak = ParaGrid.ThetaList(kTheta);
            rGrid(:, kTheta) = (  (b(kTheta)) ./ sqrt(1-(e(kTheta).*cos(ParaGrid.PhiList)).^2) );
            [Y, Z] = pol2cart( ParaGrid.PhiList, rGrid(:, kTheta) );
%             X = 0*rGrid(:, kTheta) + rGrid(end, kTheta) ./ tan(Thetak);
            X = rGrid(:, kTheta) ./ tan(Thetak);
            rGrid(:, kTheta) = sqrt(X.^2 + Y.^2 + Z.^2);
    end
    
    rGrid (:, 1) = norm(ParaSystem.Nose.KSMU_Rp .* (ParaSystem.Rp ./ (ParaSystem.r0)));
    rSol = rGrid;
    
    
%% TEST: NMM
%}
    
    
%         rGrid = rGrid_ini;

    

%}


%{
% Optimisation Grid


    N = length(ParaGrid.PhiList);
    
    I1 = (1:length(ParaGrid.PhiList)*length(ParaGrid.ThetaList)).';
    J1 = I1;
    V1 = 0*I1 + 1;
    
    I2 = I1(:, 2:end);
    J2 = I2-1;
    V2 = 0*I2 + 1;
    
    I3 = I1(:, 1:end-1);
    J3 = I3+1;
    V3 = 0*I3+1;
    
    I4 = I1(:, N+1:end);
    J4 = I4-N;
    V4 = 0*I4+1;
    
    I1Flipped = flipud(I1);
    I1Flipped = I1Flipped(N+1:end);
    I5 = flipud(I1Flipped);
    J5 = I5+N;
    V5 = 0*I5+1;
    
    I = vertcat(I1, I2, I3, I4, I5);
    J = vertcat(J1, J2, J3, J4, J5);
    V = vertcat(V1, V2, V3, V4, V5);

    Jstr = sparse(I, J, V);

%}


%{

    %% TEST: Comparing Jacobian
    
%        [PBGrid, JacobianMatrix] = PB_Grid(rGrid, ParaGrid, ParaSystem);
%        figure;
%        ErrorJac = abs(JacobianMatrix - grad_fd)
%        pcolor(ErrorJac);
%        colorbar
%        caxis([0 1.e-1])

%}