function SystemParameters = SystemParameters(Planet)

    if Planet == "Saturn"

        % Planet and Pressure Balance at SATURN
    	SystemParameters.beta = 2;             % Plasma Beta (particle pressure / magnetic pressure)
    	SystemParameters.Bp = 20000*10^(-9);   % Equatorial field, T
        SystemParameters.Rp = 60280*10^3;      % Planet radius, m
        mu0 = 4*pi*10^(-7);
        Pswnpa = 0.02;
        SystemParameters.b0 = sqrt(2*mu0*Pswnpa*10^(-9));                                                      % Field scale
        SystemParameters.r0 = (2*SystemParameters.Bp*SystemParameters.Rp^3/SystemParameters.b0)^(1/3);         % Distance scale
        SystemParameters.M = SystemParameters.Bp*SystemParameters.Rp^3 / (SystemParameters.b0*SystemParameters.r0^3);

        % CAN Disk Parameters (Bunce2007)
        SystemParameters.mu0I = ( -60.4*10^(-9) /SystemParameters.b0 ) ;      % (mu0 I0 = 60.4 nT)
        SystemParameters.a = 8;      % Inner radius of the current distribution (in Rp)
        SystemParameters.b = 15.5;   % Outer radius of the current distribution (in Rp)
        SystemParameters.D = 3;      % Disk Semi-thickness (in Rp)
        SystemParameters.CAN_DiskParameters = [SystemParameters.mu0I SystemParameters.a SystemParameters.b SystemParameters.D];

    elseif Planet == "Jupiter"

        % Planet and Pressure Balance at JUPITER
        SystemParameters.beta = 15;              % Plasma Beta (particle pressure / magnetic pressure)
        SystemParameters.Bp = 776.6*10^(-6);     % Equatorial field, T
        SystemParameters.Rp = 69911*10^3;        % Planet radius, m
        nJ = 0.4;                                % particle/cm3
        vJ = 400;                                % km/s
        mp = 1.6726*10^(-6);                     % proton mass
        PswnpaJ = mp*nJ * vJ^2;
        mu0 = 4*pi*10^(-7);
        SystemParameters.b0 = sqrt(2*mu0*PswnpaJ*10^(-9));                                                     % Field scale
        SystemParameters.r0 = (2*SystemParameters.Bp*SystemParameters.Rp^3/SystemParameters.b0)^(1/3);         % Distance scale
        SystemParameters.M = SystemParameters.Bp*SystemParameters.Rp^3 / (SystemParameters.b0*SystemParameters.r0^3);

        % CAN Disk Parameters JUPITER (Connerney1981)
        SystemParameters.mu0I = ( -(17.1*10^6/SystemParameters.Rp*mu0) /SystemParameters.b0 ) ;      % (mu0 I0 = 52.1 nT)
        SystemParameters.a = 5;      % Inner radius of the current distribution (in Rp)
        SystemParameters.b = 30;   % Outer radius of the current distribution (in Rp)
        SystemParameters.D = 2.5;      % Disk Semi-thickness (in Rp)
        SystemParameters.CAN_DiskParameters = [SystemParameters.mu0I SystemParameters.a SystemParameters.b SystemParameters.D];

    end

end
