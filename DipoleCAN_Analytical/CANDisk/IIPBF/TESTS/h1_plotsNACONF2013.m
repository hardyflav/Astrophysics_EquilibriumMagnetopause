% h1 Plots
% generate plots of h1 for Case 35 and Case 23 
% (see slide 15 of the presentation)
% you must modify IIPBF to output the h1 function and 
% first two zeros since they are not normal outputs

% Reference:
% Ratnanather, J. T., Kim, J. H., Zhang, S., Davis, A. M. J., and Lucas, 
% S. K. 2013. Algorithm XXX: IIPBF, a MATLAB toolbox for infinite integral 
% of products of two Bessel functions. To appear in ACM Transactions on 
% Mathematical Software. 

%% Case 35
% high error case

addpath('..')

abserr = 1e-10; relerr = 1e-11;

rho = 10;
tau = 1001;
c = 2.5;
a = 4.2;
b = 0.3;
fx = @(x) x.^(-c);
type = 'JY';

[result,err,neval]=IIPBF(fx,rho,tau,a,b,abserr,relerr,type);

answer = (2./pi) .* sin(pi.*(a-b-c)./2) ...
    .* ((rho.^a.*gamma(.5-.5.*c+.5.*b+.5.*a).*gamma(.5-.5.*c-.5.*b+.5.*a)) ...
	./ (2.^(c+1).*tau.^(-c+a+1).*gamma(a+1))) ...
	.* hypergeom([.5-.5.*c+.5.*b+.5.*a, .5-.5.*c-.5.*b+.5.*a], a+1, ((rho.^2)./(tau.^2)));

error=abs(result-answer);

load h1plots.mat
xmin = 0; xmax = 2.5*zero2 - zero1;
ymin = -1e7; ymax = 1e7;
figure
fplot(h1, [xmin xmax ymin ymax])
set(findobj('Type','line'),'Linewidth',2)
hold on
scatter([zero1 zero2],[h1(zero1) h1(zero2)],'k','fill')
set(gca,'Fontsize',16,'Xtick',[0 0.005 0.01 0.015],'Ytick',[-1e7 -5e6 0 5e6 1e7])
xlabel('x');ylabel('h1')
text(zero1,h1(zero1),'Zero 1','VerticalAlign','bottom','HorizontalAlign','right','Fontsize',16)
text(zero2,h1(zero2),'Zero 2','VerticalAlign','bottom','HorizontalAlign','left','Fontsize',16)
hold off

%%
% low error case

abserr = 1e-10; relerr = 1e-11;

rho = 10;
tau = 1001;
c = 1;
a = 2.5;
b = 1.5;
fx = @(x) x.^(-c);
type = 'JY';

[result,err,neval]=IIPBF(fx,rho,tau,a,b,abserr,relerr,type);

answer = (2./pi) .* sin(pi.*(a-b-c)./2) ...
    .* ((rho.^a.*gamma(.5-.5.*c+.5.*b+.5.*a).*gamma(.5-.5.*c-.5.*b+.5.*a)) ...
	./ (2.^(c+1).*tau.^(-c+a+1).*gamma(a+1))) ...
	.* hypergeom([.5-.5.*c+.5.*b+.5.*a, .5-.5.*c-.5.*b+.5.*a], a+1, ((rho.^2)./(tau.^2)));

error=abs(result-answer);

load h1plots.mat
xmin = 0; xmax = 2.5*zero2 - zero1;
ymin = -1e3; ymax = 1e3;
figure
fplot(h1, [xmin xmax ymin ymax])
set(findobj('Type','line'),'Linewidth',2)
hold on
scatter([zero1 zero2],[h1(zero1) h1(zero2)],'k','fill')
set(gca,'Fontsize',16,'Xtick',[0 0.005 0.01 0.015],'Ytick',[-1000 -500 0 500 1000])
xlabel('x');ylabel('h1')
text(zero1,h1(zero1),'Zero 1','VerticalAlign','bottom','HorizontalAlign','right','Fontsize',16)
text(zero2,h1(zero2),'Zero 2','VerticalAlign','bottom','HorizontalAlign','left','Fontsize',16)
hold off

%% Case 23
% high error case

abserr = 1e-10; relerr = 1e-11;

rho = 0.001;
tau = 0.0011;
c = 2;
a = 3;
b = 5;
fx = @(x)power(x,b-a+1)./(c.^2+x.^2);
type = 'JJ';

[result,err,neval]=IIPBF(fx,rho,tau,a,b,abserr,relerr,type);

answer = c^(b-a)*besseli(a,rho*c)*besselk(b,tau*c);

error=abs(result-answer);

load h1plots.mat
xmin = 0; xmax = 2.5*zero2 - zero1;
ymin = -0.25; ymax = 0.25;
figure
fplot(h1, [xmin xmax ymin ymax])
set(findobj('Type','line'),'Linewidth',2)
hold on
scatter([zero1 zero2],[h1(zero1) h1(zero2)],'k','fill')
set(gca,'Fontsize',16,'Xtick',[0 5000 10000 15000 20000],'Ytick',[-0.2 -0.1 0 0.1 0.2])
xlabel('x');ylabel('h1')
text(zero1,h1(zero1),'Zero 1','VerticalAlign','bottom','HorizontalAlign','right','Fontsize',16)
text(zero2,h1(zero2),'Zero 2','VerticalAlign','bottom','HorizontalAlign','left','Fontsize',16)
hold off

%% 
% low error case

abserr = 1e-10; relerr = 1e-11;

rho = 0.001;
tau = 0.0011;
c = 2;
a = 1;
b = 1;
fx = @(x)power(x,b-a+1)./(c.^2+x.^2);
type = 'JJ';

[result,err,neval]=IIPBF(fx,rho,tau,a,b,abserr,relerr,type);

answer = c^(b-a)*besseli(a,rho*c)*besselk(b,tau*c);

error=abs(result-answer);

load h1plots.mat
xmin = 0; xmax = 2.5*zero2 - zero1;
ymin = -0.25; ymax = 0.25;
figure
fplot(h1, [xmin xmax ymin ymax])
set(findobj('Type','line'),'Linewidth',2)
hold on
scatter([zero1 zero2],[h1(zero1) h1(zero2)],'k','fill')
set(gca,'Fontsize',16,'Xtick',[0 2500 5000 7500 10000],'Ytick',[-0.2 -0.1 0 0.1 0.2])
xlabel('x');ylabel('h1')
text(zero1,h1(zero1),'Zero 1','VerticalAlign','top','HorizontalAlign','right','Fontsize',16)
text(zero2,h1(zero2),'Zero 2','VerticalAlign','top','HorizontalAlign','left','Fontsize',16)
hold off

