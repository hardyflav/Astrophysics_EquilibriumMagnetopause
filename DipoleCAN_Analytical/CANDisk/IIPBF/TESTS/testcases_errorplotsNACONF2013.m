% run case 6 with complex-valued c, cases 33-35 and plot Actual error and
% Estimated error vs. # of evaluations to evaluate performance and accuracy
% of IIPBF.
% 
%   Reference: Ratnanather, J. T., Kim, J. H., Zhang, S., Davis, A. M. J., 
%              and Lucas, S. K. 2012. IIPBF: a MATLAB toolbox for infinite 
%              integrals of product of Bessel functions. Pending review
%              from ACM Transactions on Mathematical Software
% 
% June 4, 2013, Geoffrey Gunter

addpath('..')

clear all
close all

% Case 6
% x * K_0(x*(10+2i)) * J_0(x) * J_0(x)
a = 10; B = 2;
c1 = a+1i*B;
rho = 0.1; tau = 100;
[sol06_1,ans06_1,err06_1,n06_1,~] = testcasesNACONF2013(6,c1,rho,tau);

abs06_1 = abs(sol06_1-ans06_1);

% figure
% fplot(@(x) x .* besselk(0,x*c1) .* besselj(0,rho.*x) .* besselj(0,tau.*x), [1 10])

fprintf('\nCase 6_a\n')
fprintf('c = %g + %gi\n',a,B)
fprintf('rho = %g, tau = %g\n',rho,tau)

% Case 6
% x * K_0(x*(2+2i)) * J_0(x) * J_0(x)
a = 2; B = 2;
c2 = a+1i*B;
rho = 0.1; tau = 100;
[sol06_2,ans06_2,err06_2,n06_2,~] = testcasesNACONF2013(6,c2,rho,tau);

abs06_2 = abs(sol06_2-ans06_2);

% figure
% fplot(@(x) x .* besselk(0,x*c2) .* besselj(0,rho.*x) .* besselj(0,tau.*x), [1 10])

fprintf('\nCase 6_b\n')
fprintf('c = %g + %gi\n',a,B)
fprintf('rho = %g, tau = %g\n',rho,tau)

% Case 6
% x * K_0(x*(2+10i)) * J_0(x) * J_0(x)
a = 2; B = 10;
c3 = a+1i*B;
rho = 0.1; tau = 100;
[sol06_3,ans06_3,err06_3,n06_3,~] = testcasesNACONF2013(6,c3,rho,tau);

abs06_3 = abs(sol06_3-ans06_3);

% figure
% fplot(@(x) x .* besselk(0,x*c3) .* besselj(0,rho.*x) .* besselj(0,tau.*x), [1 10])

fprintf('\nCase 6_c\n')
fprintf('c = %g + %gi\n',a,B)
fprintf('rho = %g, tau = %g\n',rho,tau)

% Case 33
% x/(x^2+2^2) * J_1(2x) * J_1(x)
c = 2;
a = 1;
rho = 2; tau = 1;
[sol33,ans33,err33,n33,~] = testcasesNACONF2013(33,c,a,rho,tau);

abs33 = abs(sol33-ans33);

fprintf('\nCase 33\n')
fprintf('c = %g\n',c)
fprintf('a = %g, b = a\n',a)
fprintf('rho = %g, tau = %g\n',rho,tau)

% Case 34
c = 2;
a = 1; b = 2;
rho = 2; tau = 1;
[sol34,ans34,err34,n34,~] = testcasesNACONF2013(34,c,a,b,rho,tau);

abs34 = abs(sol34-ans34);

fprintf('\nCase 34\n')
fprintf('c = %g\n',c)
fprintf('a = %g, b = %g\n',a,b)
fprintf('rho = %g, tau = %g\n',rho,tau)

% Case 35
c = 1;
a = 2; b = 1;
rho = 1; tau = 2;
[sol35,ans35,err35,n35,~] = testcasesNACONF2013(35,c,a,b,rho,tau);

abs35 = abs(sol35-ans35);

fprintf('\nCase 35\n')
fprintf('c = %g\n',c)
fprintf('a = %g, b = %g\n',a,b)
fprintf('rho = %g, tau = %g\n',rho,tau)


figure

semilogy (n06_1, abs06_1, 'k-.','Linewidth',3)
hold on
semilogy (n06_2, abs06_2, 'ys-','Linewidth',3)
semilogy (n06_3, abs06_3, 'r*-')
semilogy (n33, abs33, 'gx-','Linewidth',2);
semilogy (n34, abs34, 'kd-.','Markersize',12)
semilogy (n35, abs35,  'bo:','Linewidth',2)

set(gca,'Fontsize',14,'Fontname','Helvetica','Linewidth',3,'Xlim',[20 80],...
    'Xtick',[0 20 40 60 80 100],'Ylim',[0 1e-8],'Ytick',[1e-22 1e-20 1e-18 1e-16 1e-14 1e-12 1e-10 1e-8])
xlabel('Number of Evaluations','Fontsize',14,'Fontname','Helvetica')
ylabel('Actual Error','Fontsize',14,'Fontname','Helvetica')
legend('6_a','6_b','6_c','33','34','35','Location','NorthEast')
print -djpeg Act_error.jpg

hold off

figure

semilogy (n06_1, err06_1, 'k-.','Linewidth',3)
hold on
semilogy (n06_2, err06_2, 'ys-','Linewidth',3)
semilogy (n06_3, err06_3, 'r*-')
semilogy (n33, err33, 'gx-','Linewidth',2);
semilogy (n34, err34, 'kd-.','Markersize',12)
semilogy (n35, err35,  'bo:','Linewidth',2)

set(gca,'Fontsize',14,'Fontname','Helvetica','Linewidth',3,'Xlim',[20 80],...
    'Xtick',[0 20 40 60 80 100],'Ylim',[0 1e-8],'Ytick',[1e-16 1e-14 1e-12 1e-10 1e-8 1e-4])
xlabel('Number of Evaluations','Fontsize',14,'Fontname','Helvetica')
ylabel('Estimated Error','Fontsize',14,'Fontname','Helvetica')
legend('6_a','6_b','6_c','33','34','35','Location','NorthEast')
print -djpeg Est_error.jpg

hold off


