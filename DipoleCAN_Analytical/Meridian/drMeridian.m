function drMeridian = drMeridian(theta, r)

    Bp = 20000*10^(-9);     % Equatorial field, T
    Rp = 60280*10^3;        % Planet radius, m
    
    mu_0 = 4*pi*10^(-7);
    Pswnpa = 0.02;
    
    b0 = sqrt(2*mu_0*Pswnpa*10^(-9));   % Field scale
    r0 = (2*Bp*Rp^3/b0)^(1/3);          % Distance scale
    M = (Bp/b0)*(Rp/r0)^3;

drMeridian = (4.*M+r.^3).^(-1).*(2.*M.*r.*cos(theta)+(-1).*r.^4.*cos(theta)).* ...
  csc(theta);

end