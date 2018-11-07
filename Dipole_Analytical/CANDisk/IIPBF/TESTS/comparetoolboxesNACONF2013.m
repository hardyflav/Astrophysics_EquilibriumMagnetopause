% compare error and computation time for toolboxes for case 33 and 34, 
% IIPBF and BESSELINT
% 
%   Reference: Ratnanather, J. T., Kim, J. H., Zhang, S., Davis, A. M. J., 
%              and Lucas, S. K. 2012. IIPBF: a MATLAB toolbox for infinite 
%              integrals of product of Bessel functions. Pending review
%              from ACM Transactions on Mathematical Software
% 
% June 6, 2013, Geoffrey Gunter

addpath('..')
addpath('./BESSELINT')

%% Case 33
% f(x) = x/(x^2+c^2)
% B = J_a(rho*x)*J_a(tau*x)
% [Re(c) > 0, Re(a) > -1]

clear all; close all;
format long e;

abserr = 1.e-14;
relerr = 1.e-13;

c = 2;
a = 1; b = a;
rho = 2; tau = 1;
fx = @(x) x./(x.^2+c.^2);
type = 'JJ';

tstart = tic;
[f33,err33] = besselintr([rho tau],[a b],1,c,relerr,abserr);
telapsed33_1 = toc(tstart);

tstart = tic;
[res33,rel33,evals33]=IIPBF(fx,rho,tau,a,b,abserr,relerr,type);
telapsed33_2 = toc(tstart);

err33_1 = err33(1)/(10.^floor(log10(err33(1))));
err33_2 = floor(log10(err33(1)));
rel33_1 = rel33/(10.^floor(log10(rel33)));
rel33_2 = floor(log10(rel33));
diff33 = abs(res33-f33);
diff33_1 = diff33/(10.^floor(log10(diff33)));
diff33_2 = floor(log10(diff33));
if (diff33==0)
    diff33_1 = 0;
    diff33_2 = 0;
end

fprintf('33 & %1.4f & %1.4f & %1.2f %s & %1.2f %s & %1.2f %s \n', telapsed33_1, telapsed33_2, err33_1, ['x$10^{', num2str(err33_2),'}$'], rel33_1, ['x$10^{', num2str(rel33_2),'}$'], diff33_1, ['x$10^{', num2str(diff33_2),'}$ \\']);

%% Case 34
% f(x) = x^-c
% B = J_a(rho*x)*J_b(tau*x)
% [Re(a+b-c+1) > 0, Re(c) > -1, tau < rho]

clear all; close all;
format long e;

abserr = 1.e-14;
relerr = 1.e-13;

c = 3;
a = 1; b = 2;
rho = 2; tau = 1;
fx = @(x) x.^(-c);
type = 'JJ';

tstart = tic;
[f34,err34] = besselint([rho tau],[a b],-c,relerr,abserr);
telapsed34_1 = toc(tstart);

tstart = tic;
[res34,rel34,evals34]=IIPBF(fx,rho,tau,a,b,abserr,relerr,type);
telapsed34_2 = toc(tstart);

err34_1 = err34(1)/(10.^floor(log10(err34(1))));
err34_2 = floor(log10(err34(1)));
rel34_1 = rel34/(10.^floor(log10(rel34)));
rel34_2 = floor(log10(rel34));
diff34 = abs(res34-f34);
diff34_1 = diff34/(10.^floor(log10(diff34)));
diff34_2 = floor(log10(diff34));
if (diff34==0)
    diff34_1 = 0;
    diff34_2 = 0;
end

fprintf('34 & %1.4f & %1.4f & %1.2f %s & %1.2f %s & %1.2f %s \n', telapsed34_1, telapsed34_2, err34_1, ['x$10^{', num2str(err34_2),'}$'], rel34_1, ['x$10^{', num2str(rel34_2),'}$'], diff34_1, ['x$10^{', num2str(diff34_2),'}$ \\']);
