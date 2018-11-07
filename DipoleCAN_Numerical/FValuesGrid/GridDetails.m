function [NbPointsTotal, NbPointsInside, NThetaInside, NPhiInside, ThetaSpanVector, PhiSpanVector] = GridDetails(ThetaMaxDeg, PhiMaxDeg, DeltaThetaDeg, DeltaPhiDeg)

    NThetaTotal = ThetaMaxDeg/DeltaThetaDeg+1;
    NPhiTotal = PhiMaxDeg/DeltaPhiDeg+1;
    NbPointsTotal = NThetaTotal*NPhiTotal;

    NThetaInside = NThetaTotal-2;
    NPhiInside = NPhiTotal-2;
    NbPointsInside = NThetaInside*NPhiInside;

    ThetaSpanVectorTemp = repmat(DeltaThetaDeg:DeltaThetaDeg:ThetaMaxDeg-DeltaThetaDeg, [NPhiInside 1]);
    ThetaSpanVector = ThetaSpanVectorTemp(:)*pi/180;

    PhiSpanVectorTemp = repmat(DeltaPhiDeg:DeltaPhiDeg:PhiMaxDeg-DeltaPhiDeg, [1, NThetaInside]).';
    PhiSpanVector = PhiSpanVectorTemp(:)*pi/180;

end
