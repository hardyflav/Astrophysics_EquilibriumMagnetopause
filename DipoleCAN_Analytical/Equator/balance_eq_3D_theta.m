function balance_eq_3D_theta = balance_eq_3D_theta (theta, phi, r, dr)

    Bp = 20000*10^(-9);     % Equatorial field, T
    Rp = 60280*10^3;        % Planet radius, m
    
    mu_0 = 4*pi*10^(-7);
    Pswnpa = 0.02;
    
    b0 = sqrt(2*mu_0*Pswnpa*10^(-9));   % Field scale
    r0 = (2*Bp*Rp^3/b0)^(1/3);          % Distance scale
    M = (Bp/b0)*(Rp/r0)^3;
    
    beta = 1;
    
n_vec = [ 1 ; (-1/r)*dr ; 0 ];
B = (M/r^3) .* [ -2*sin(theta)*sin(phi) ; cos(theta)*sin(phi) ; cos(phi) ];
v = [ -cos(theta) ; sin(theta) ; 0 ];

balance_eq_3D_theta = norm(cross(n_vec, B)) + (1/2)*dot(n_vec,v);

end
