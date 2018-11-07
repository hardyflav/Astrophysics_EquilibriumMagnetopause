function [NbPointsTotal, NbPointsInside, NThetaInside, NPhiInside, ThetaSpanVector, PhiSpanVector] = GridDetailsBottom(ThetaLeftBoundaryDeg, ThetaMaxDeg, PhiMaxDeg, DeltaThetaDeg, DeltaPhiDeg)

ThetaMinDeg = ThetaLeftBoundaryDeg + DeltaThetaDeg;

NThetaTotal = length(ThetaLeftBoundaryDeg:DeltaThetaDeg:ThetaMaxDeg);
NPhiTotal = length(0:DeltaPhiDeg:PhiMaxDeg);
NbPointsTotal = NThetaTotal*NPhiTotal;

NThetaInside = NThetaTotal-2;
NPhiInside = NPhiTotal-2;
NbPointsInside = NThetaInside*NPhiInside;

ThetaSpanVectorTemp = repmat(ThetaMinDeg:DeltaThetaDeg:ThetaMaxDeg-DeltaThetaDeg, [NPhiInside 1]);
ThetaSpanVector = ThetaSpanVectorTemp(:)*pi/180;

PhiSpanVectorTemp = repmat(DeltaPhiDeg:DeltaPhiDeg:PhiMaxDeg-DeltaPhiDeg, [1, NThetaInside]).';
PhiSpanVector = PhiSpanVectorTemp(:)*pi/180;


end
