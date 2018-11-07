% parameter tests for Case 35, 
% generate error tables for slide 14 of presentation

% Reference:
% Ratnanather, J. T., Kim, J. H., Zhang, S., Davis, A. M. J., and Lucas, 
% S. K. 2013. Algorithm XXX: IIPBF, a MATLAB toolbox for infinite integral 
% of products of two Bessel functions. To appear in ACM Transactions on 
% Mathematical Software. 

addpath('..')
warning off

%%
% a = 2.5, b = 1.5, c = 1
rho_vec = [0.001 0.01 0.1 1 10 100 1000];
tau_vec = [0.0011 0.011 0.11 1.1 11 101 1001];
% tau > rho for case 35
abserr = 1e-10; relerr = 1e-11;

c = 1;
a = 2.5; b = 1.5;
fx = @(x) x.^(-c);
type = 'JY';

for i=1:length(rho_vec)
    for j=1:length(tau_vec)
        if rho_vec(i)<tau_vec(j)
            rho=rho_vec(i); tau=tau_vec(j);
            % fprintf('rho = %g, tau = %g\n',rho,tau)
            [results(i,j),rel_err(i,j),neval(i,j)]=IIPBF(fx,rho,tau,a,b,abserr,relerr,type);
            
            answers(i,j) = (2./pi) .* sin(pi.*(a-b-c)./2) ...
                .* ((rho.^a.*gamma(.5-.5.*c+.5.*b+.5.*a).*gamma(.5-.5.*c-.5.*b+.5.*a)) ...
                ./ (2.^(c+1).*tau.^(-c+a+1).*gamma(a+1))) ...
                .* hypergeom([.5-.5.*c+.5.*b+.5.*a, .5-.5.*c-.5.*b+.5.*a], a+1, ((rho.^2)./(tau.^2)));
            
            abs_err(i,j) = abs(results(i,j) - answers(i,j));
            
            act_err_1(i,j) = abs_err(i,j)/(10.^floor(log10(abs_err(i,j))));
            act_err_2(i,j) = floor(log10(abs_err(i,j)));
            if (abs_err(i,j)==0)
                act_err_1(i,j) = 0;
                act_err_2(i,j) = 0;
            end
        end
    end
end

fprintf('\n\\multicolumn{9}{c}{$a=%g, b=%g, c=%g$} \\\\\n', a, b, c);
fprintf('       &        &               \\multicolumn{7}{c}{$\\rho$} \\\\\n');
fprintf('       &        & %s & %s & %s & %s & %s & %s & %s \\\\ \\cline{3-9}\n', num2str(rho_vec(1)), num2str(rho_vec(2)), num2str(rho_vec(3)), num2str(rho_vec(4)), num2str(rho_vec(5)), num2str(rho_vec(6)), num2str(rho_vec(7)));
fprintf('       & %s & %1.2f %s & - & - & - & - & - & - \\\\\n', num2str(tau_vec(1)), act_err_1(1,1), ['x$10^{' num2str(act_err_2(1,1))  '}$']);
fprintf('       & %s & %1.2f %s & %1.2f %s & - & - & - & - & - \\\\\n', num2str(tau_vec(2)), act_err_1(1,2), ['x$10^{' num2str(act_err_2(1,2))  '}$'],  act_err_1(2,2), ['x$10^{' num2str(act_err_2(2,2))  '}$']);
fprintf('       & %s & %1.2f %s & %1.2f %s & %1.2f %s & - & - & - & - \\\\\n', num2str(tau_vec(3)), act_err_1(1,3), ['x$10^{' num2str(act_err_2(1,3))  '}$'],  act_err_1(2,3), ['x$10^{' num2str(act_err_2(2,3))  '}$'], act_err_1(3,3), ['x$10^{' num2str(act_err_2(3,3))  '}$']);
fprintf('$\\tau$ & %s & %1.2f %s & %1.2f %s & %1.2f %s & %1.2f %s & - & - & - \\\\\n', num2str(tau_vec(4)), act_err_1(1,4), ['x$10^{' num2str(act_err_2(1,4))  '}$'],  act_err_1(2,4), ['x$10^{' num2str(act_err_2(2,4))  '}$'], act_err_1(3,4), ['x$10^{' num2str(act_err_2(3,4))  '}$'], act_err_1(4,4), ['x$10^{' num2str(act_err_2(4,4))  '}$']);
fprintf('       & %s & %1.2f %s & %1.2f %s & %1.2f %s & %1.2f %s & %1.2f %s & - & - \\\\\n', num2str(tau_vec(5)), act_err_1(1,5), ['x$10^{' num2str(act_err_2(1,5))  '}$'],  act_err_1(2,5), ['x$10^{' num2str(act_err_2(2,5))  '}$'], act_err_1(3,5), ['x$10^{' num2str(act_err_2(3,5))  '}$'], act_err_1(4,5), ['x$10^{' num2str(act_err_2(4,5))  '}$'], act_err_1(5,5), ['x$10^{' num2str(act_err_2(5,5))  '}$']);
fprintf('       & %s & %1.2f %s & %1.2f %s & %1.2f %s & %1.2f %s & %1.2f %s & %1.2f %s & - \\\\\n', num2str(tau_vec(6)), act_err_1(1,6), ['x$10^{' num2str(act_err_2(1,6))  '}$'],  act_err_1(2,6), ['x$10^{' num2str(act_err_2(2,6))  '}$'], act_err_1(3,6), ['x$10^{' num2str(act_err_2(3,6))  '}$'], act_err_1(4,6), ['x$10^{' num2str(act_err_2(4,6))  '}$'], act_err_1(5,6), ['x$10^{' num2str(act_err_2(5,6))  '}$'], act_err_1(6,6), ['x$10^{' num2str(act_err_2(6,6))  '}$']);
fprintf('       & %s & %1.2f %s & %1.2f %s & %1.2f %s & %1.2f %s & %1.2f %s & %1.2f %s & %1.2f %s \\\\ \\hline\n', num2str(tau_vec(7)), act_err_1(1,7), ['x$10^{' num2str(act_err_2(1,7))  '}$'],  act_err_1(2,7), ['x$10^{' num2str(act_err_2(2,7))  '}$'], act_err_1(3,7), ['x$10^{' num2str(act_err_2(3,7))  '}$'], act_err_1(4,7), ['x$10^{' num2str(act_err_2(4,7))  '}$'], act_err_1(5,7), ['x$10^{' num2str(act_err_2(5,7))  '}$'], act_err_1(6,7), ['x$10^{' num2str(act_err_2(6,7))  '}$'], act_err_1(7,7), ['x$10^{' num2str(act_err_2(7,7))  '}$']);

%%
% a = 4.2, b = 0.3, c = 2.5
rho_vec = [0.001 0.01 0.1 1 10 100 1000];
tau_vec = [0.0011 0.011 0.11 1.1 11 101 1001];
% tau > rho for case 35
abserr = 1e-10; relerr = 1e-11;

c = 2.5;
a = 4.2; b = 0.3;
fx = @(x) x.^(-c);
type = 'JY';

for i=1:length(rho_vec)
    for j=1:length(tau_vec)
        if rho_vec(i)<tau_vec(j)
            rho=rho_vec(i); tau=tau_vec(j);
            % fprintf('rho = %g, tau = %g\n',rho,tau)
            [results(i,j),rel_err(i,j),neval(i,j)]=IIPBF(fx,rho,tau,a,b,abserr,relerr,type);
            
            answers(i,j) = (2./pi) .* sin(pi.*(a-b-c)./2) ...
                .* ((rho.^a.*gamma(.5-.5.*c+.5.*b+.5.*a).*gamma(.5-.5.*c-.5.*b+.5.*a)) ...
                ./ (2.^(c+1).*tau.^(-c+a+1).*gamma(a+1))) ...
                .* hypergeom([.5-.5.*c+.5.*b+.5.*a, .5-.5.*c-.5.*b+.5.*a], a+1, ((rho.^2)./(tau.^2)));
            
            abs_err(i,j) = abs(results(i,j) - answers(i,j));
            
            act_err_1(i,j) = abs_err(i,j)/(10.^floor(log10(abs_err(i,j))));
            act_err_2(i,j) = floor(log10(abs_err(i,j)));
            if (abs_err(i,j)==0)
                act_err_1(i,j) = 0;
                act_err_2(i,j) = 0;
            end
        end
    end
end

fprintf('\n\\multicolumn{9}{c}{$a=%g, b=%g, c=%g$} \\\\\n', a, b, c);
fprintf('       &        &               \\multicolumn{7}{c}{$\\rho$} \\\\\n');
fprintf('       &        & %s & %s & %s & %s & %s & %s & %s \\\\ \\cline{3-9}\n', num2str(rho_vec(1)), num2str(rho_vec(2)), num2str(rho_vec(3)), num2str(rho_vec(4)), num2str(rho_vec(5)), num2str(rho_vec(6)), num2str(rho_vec(7)));
fprintf('       & %s & %1.2f %s & - & - & - & - & - & - \\\\\n', num2str(tau_vec(1)), act_err_1(1,1), ['x$10^{' num2str(act_err_2(1,1))  '}$']);
fprintf('       & %s & %1.2f %s & %1.2f %s & - & - & - & - & - \\\\\n', num2str(tau_vec(2)), act_err_1(1,2), ['x$10^{' num2str(act_err_2(1,2))  '}$'],  act_err_1(2,2), ['x$10^{' num2str(act_err_2(2,2))  '}$']);
fprintf('       & %s & %1.2f %s & %1.2f %s & %1.2f %s & - & - & - & - \\\\\n', num2str(tau_vec(3)), act_err_1(1,3), ['x$10^{' num2str(act_err_2(1,3))  '}$'],  act_err_1(2,3), ['x$10^{' num2str(act_err_2(2,3))  '}$'], act_err_1(3,3), ['x$10^{' num2str(act_err_2(3,3))  '}$']);
fprintf('$\\tau$ & %s & %1.2f %s & %1.2f %s & %1.2f %s & %1.2f %s & - & - & - \\\\\n', num2str(tau_vec(4)), act_err_1(1,4), ['x$10^{' num2str(act_err_2(1,4))  '}$'],  act_err_1(2,4), ['x$10^{' num2str(act_err_2(2,4))  '}$'], act_err_1(3,4), ['x$10^{' num2str(act_err_2(3,4))  '}$'], act_err_1(4,4), ['x$10^{' num2str(act_err_2(4,4))  '}$']);
fprintf('       & %s & %1.2f %s & %1.2f %s & %1.2f %s & %1.2f %s & %1.2f %s & - & - \\\\\n', num2str(tau_vec(5)), act_err_1(1,5), ['x$10^{' num2str(act_err_2(1,5))  '}$'],  act_err_1(2,5), ['x$10^{' num2str(act_err_2(2,5))  '}$'], act_err_1(3,5), ['x$10^{' num2str(act_err_2(3,5))  '}$'], act_err_1(4,5), ['x$10^{' num2str(act_err_2(4,5))  '}$'], act_err_1(5,5), ['x$10^{' num2str(act_err_2(5,5))  '}$']);
fprintf('       & %s & %1.2f %s & %1.2f %s & %1.2f %s & %1.2f %s & %1.2f %s & %1.2f %s & - \\\\\n', num2str(tau_vec(6)), act_err_1(1,6), ['x$10^{' num2str(act_err_2(1,6))  '}$'],  act_err_1(2,6), ['x$10^{' num2str(act_err_2(2,6))  '}$'], act_err_1(3,6), ['x$10^{' num2str(act_err_2(3,6))  '}$'], act_err_1(4,6), ['x$10^{' num2str(act_err_2(4,6))  '}$'], act_err_1(5,6), ['x$10^{' num2str(act_err_2(5,6))  '}$'], act_err_1(6,6), ['x$10^{' num2str(act_err_2(6,6))  '}$']);
fprintf('       & %s & %1.2f %s & %1.2f %s & %1.2f %s & %1.2f %s & %1.2f %s & %1.2f %s & %1.2f %s \\\\ \\hline\n', num2str(tau_vec(7)), act_err_1(1,7), ['x$10^{' num2str(act_err_2(1,7))  '}$'],  act_err_1(2,7), ['x$10^{' num2str(act_err_2(2,7))  '}$'], act_err_1(3,7), ['x$10^{' num2str(act_err_2(3,7))  '}$'], act_err_1(4,7), ['x$10^{' num2str(act_err_2(4,7))  '}$'], act_err_1(5,7), ['x$10^{' num2str(act_err_2(5,7))  '}$'], act_err_1(6,7), ['x$10^{' num2str(act_err_2(6,7))  '}$'], act_err_1(7,7), ['x$10^{' num2str(act_err_2(7,7))  '}$']);

warning on

