
function [Btot] = MagField_2(rGrid, ParaGrid, ParaSystem_ini, Btot_km1)


    global ParaSystem
    ParaSystem = ParaSystem_ini;
    ParaGrid_Whole = ParaGrid;
    
    
%% Building Entire MP Boundary, in r0
    rWhole = vertcat(rGrid, flipud( rGrid(2:end, :)) );
    ParaGrid_Whole.PhiList = vertcat(ParaGrid.PhiList, pi+ParaGrid.PhiList(2:end, :) );
    [PhiGrid, ThetaGrid] = ndgrid(ParaGrid_Whole.PhiList, ParaGrid.ThetaList);
    ParaGrid_Whole.PhiGrid = PhiGrid;
    ParaGrid_Whole.ThetaGrid = ThetaGrid;
    

%% Dip - CAN Field
    B_DipCAN = MagField_DipCAN(rGrid, ParaGrid, ParaSystem_ini);
    B_DipCAN.X = B_DipCAN.XInterp( ParaGrid_Whole.PhiGrid, ParaGrid_Whole.ThetaGrid );
    B_DipCAN.Y = B_DipCAN.YInterp( ParaGrid_Whole.PhiGrid, ParaGrid_Whole.ThetaGrid );
    B_DipCAN.Z = B_DipCAN.ZInterp( ParaGrid_Whole.PhiGrid, ParaGrid_Whole.ThetaGrid );
    
    
%% Shielding Field
   
if ParaSystem.OnOff_ShieldingField == 0
    Btot = B_DipCAN;
else
    
    [YM, ZM, XM] = sph2cart(ParaGrid_Whole.PhiGrid, pi/2 - ParaGrid_Whole.ThetaGrid, rWhole);
    XM_List = XM(:);
    YM_List = YM(:);
    ZM_List = ZM(:);
    B_ShieldX = 0.* XM_List;
    
    if not( isequal( Btot_km1, [0; 0; 0]) )
            
        for kM=1:length(XM_List)
            
            100*(kM/length(XM_List))
            
            M = [ XM_List(kM),  YM_List(kM),  ZM_List(kM) ];
            dBpX = @(X) dBpM(rGrid, [X(:,1), X(:,2)], M, ParaGrid, ParaSystem, Btot_km1, 'X');
            dBpY = @(X) dBpM(rGrid, [X(:,1), X(:,2)], M, ParaGrid, ParaSystem, Btot_km1, 'Y');
            dBpZ = @(X) dBpM(rGrid, [X(:,1), X(:,2)], M, ParaGrid, ParaSystem, Btot_km1, 'Z');
            
            [QX,~, ~, ~] = integralN_mc( dBpX, [0, 2*pi; 0, ParaGrid.ThetaList(end)]  )   ; 
                B_ShieldX(kM) = QX;
            [QY,~, ~, ~] = integralN_mc( dBpY, [0, 2*pi; 0, ParaGrid.ThetaList(end)]  )   ; 
                B_ShieldY(kM) = QY;
            [QZ,~, ~, ~] = integralN_mc( dBpZ, [0, 2*pi; 0, ParaGrid.ThetaList(end)]  )   ; 
                B_ShieldZ(kM) = QZ;
        end
        
    else
        
        for kM=1:length(XM_List)
            
            100*(kM/length(XM_List))

            M = [ XM_List(kM),  YM_List(kM),  ZM_List(kM) ];
            dBpX = @(X) dBpM(rGrid, [X(:,1), X(:,2)], M, ParaGrid, ParaSystem, B_DipCAN, 'X');
            dBpY = @(X) dBpM(rGrid, [X(:,1), X(:,2)], M, ParaGrid, ParaSystem, B_DipCAN, 'Y');
            dBpZ = @(X) dBpM(rGrid, [X(:,1), X(:,2)], M, ParaGrid, ParaSystem, B_DipCAN, 'Z');
            
            [QX,~, ~, ~] = integralN_mc( dBpX, [0, 2*pi; 0, ParaGrid.ThetaList(end)]  )   ; 
                B_ShieldX(kM) = QX;
            [QY,~, ~, ~] = integralN_mc( dBpY, [0, 2*pi; 0, ParaGrid.ThetaList(end)]  )   ; 
                B_ShieldY(kM) = QY;
            [QZ,~, ~, ~] = integralN_mc( dBpZ, [0, 2*pi; 0, ParaGrid.ThetaList(end)]  )   ; 
                B_ShieldZ(kM) = QZ;
        end
       
        
    end
    
    BShield.X = reshape(B_ShieldX, size(rWhole));
    BShield.Y = reshape(B_ShieldY, size(rWhole));
    BShield.Z = reshape(B_ShieldZ, size(rWhole));

end
    
    
%% Packaging Output

Btot.X = B_DipCAN.X + BShield.X;
Btot.Y = B_DipCAN.Y + BShield.Y;
Btot.Z = B_DipCAN.Z + BShield.Z;

Btot.XInterp = griddedInterpolant( ParaGrid_Whole.PhiGrid, ParaGrid_Whole.ThetaGrid, Btot.X );
Btot.YInterp = griddedInterpolant( ParaGrid_Whole.PhiGrid, ParaGrid_Whole.ThetaGrid, Btot.Y );
Btot.ZInterp = griddedInterpolant( ParaGrid_Whole.PhiGrid, ParaGrid_Whole.ThetaGrid, Btot.Z );

ParaSystem.Btot_km1 = Btot;

end