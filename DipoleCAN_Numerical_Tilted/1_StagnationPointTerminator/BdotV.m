
function [BComponentV] = BdotV(CoordinateDetails, PointCoordinates, ParaSystem)


%--- Local Magnetic Field             
    B_Cartesian_KSM = BField(CoordinateDetails, PointCoordinates, ParaSystem);    

%--- Solar Wind Vector            
	v_Cartesian_KSM = [ -1 ; 0 ; 0 ];
     
    if strcmp(CoordinateDetails.Component, 'BdotV')
        
        BComponentV = dot(v_Cartesian_KSM, B_Cartesian_KSM);

    elseif strcmp(CoordinateDetails.Component, 'BcrossV')
        
        BComponentV = dot([0;0;1], B_Cartesian_KSM);
    end

end