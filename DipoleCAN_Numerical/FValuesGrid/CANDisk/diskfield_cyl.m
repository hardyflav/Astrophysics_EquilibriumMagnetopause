function Bcyl_disk=diskfield_cyl(x,xdata)
%
% provides cylindrical components of field for axisymmetric ring current
% extending from a to b, with half-thickness D, and current parameter I
%
% input: x = [I a b D]
%        xdata = [z; rho] in 2 column format (phi is irrelevant)
%
% output: Bcyl_disk, contains both components sequentially
% ie
% Bz_disk = Bcyl_disk(1:length(Bcyl_disk)/2);
% Brho_disk= Bcyl_disk(length(Bcyl_disk)/2+1:length(Bcyl_disk));
%
% Giacomo Giampieri, October 2002

current = x(1);
a = x(2);
b = x(3);
D = x(4);

z = xdata(1,:);
rho = xdata(2,:);


[Brho_disk_a,Bz_disk_a]= diskfield_inf(rho,z,a,D);
[Brho_disk_b,Bz_disk_b]= diskfield_inf(rho,z,b,D);
Brho_disk = Brho_disk_a - Brho_disk_b;
Bz_disk = Bz_disk_a - Bz_disk_b;

Bcyl_disk = current*[Bz_disk Brho_disk];
 
return;

function [Brho_disk,Bz_disk]=diskfield_inf(rho,z,a,D)
% provides cylindrical components of field for axisymmetric ring current
% extending from a to Infinity, with half-thickness D.

maxorder = 10;

tmp1 = ((z-D).^2+a^2);
tmp2 = ((z+D).^2+a^2);
tmp3 = ((z-D).^2+rho.^2);
tmp4 = ((z+D).^2+rho.^2);

for k = 0:maxorder
         
        temp1 = legendre(2*k,(z-D)./sqrt(tmp1)); 
        temp2 = legendre(2*k,(z+D)./sqrt(tmp2)); 
         
        Brho_int(:,k+1) = (((-1)^k*rho.^(2*k+1)*factorial(2*k)./(2^(2*k+2)*factorial(k)^2*(k+1))...
            .*(temp1(1,:)./tmp1.^(k+0.5) - temp2(1,:)./tmp2.^(k+0.5))))';

        if (k==0)
            
            Bz_int(:,k+1) = (1/2*log((z+D+sqrt(tmp2))./(z-D+sqrt(tmp1))))';
            
            Brho_ext(:,k+1) = (1./(2*rho).*(2*D*sign(z).*(abs(z)>=D)+2*z.*(abs(z)<D)+sqrt(tmp3)-sqrt(tmp4)))';
            
            Bz_ext(:,k+1) = (1/2*log((z+D+sqrt(tmp4))./(z-D+sqrt(tmp3))))';

        else
            
            temp3 = legendre(2*k-1,(z-D)./sqrt(tmp1)); 
            temp4 = legendre(2*k-1,(z+D)./sqrt(tmp2)); 
            temp5 = legendre(2*k-1,(z-D)./sqrt(tmp3)); 
            temp6 = legendre(2*k-1,(z+D)./sqrt(tmp4)); 

            Bz_int(:,k+1) = ((-1)^k*rho.^(2*k)*factorial(2*k-1)/(2^(2*k+1)*factorial(k)^2)...
            .*(temp3(1,:)./tmp1.^(k) - temp4(1,:)./tmp2.^(k)))';
        
            Brho_ext(:,k+1) = ((-1)^(k+1)*a^(2*k)*factorial(2*k-2)/(2^(2*k+1)*factorial(k)^2)...
            .*(temp5(2,:)./tmp3.^(k) - temp6(2,:)./tmp4.^(k)))';
       
            Bz_ext(:,k+1) = ((-1)^(k)*a^(2*k)*factorial(2*k-1)/(2^(2*k+1)*factorial(k)^2)...
            .*(temp5(1,:)./tmp3.^(k) - temp6(1,:)./tmp4.^(k)))';
        
        end
    end
    
    % now sum over k
    Brho_int_tot = sum(Brho_int,2);
    Bz_int_tot = sum(Bz_int,2);
    Brho_ext_tot = sum(Brho_ext,2);
    Bz_ext_tot = sum(Bz_ext,2);
    
    Brho_disk = Brho_int_tot'.*(rho<=a) + Brho_ext_tot'.*(rho>a);
    Bz_disk = Bz_int_tot'.*(rho<=a) + Bz_ext_tot'.*(rho>a); 
   
return;    
