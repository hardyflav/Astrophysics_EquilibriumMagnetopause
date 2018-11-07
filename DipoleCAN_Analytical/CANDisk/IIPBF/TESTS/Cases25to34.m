% error tables for Cases 25-34

% Reference:
% Ratnanather, J. T., Kim, J. H., Zhang, S., Davis, A. M. J., and Lucas, 
% S. K. 2013. Algorithm XXX: IIPBF, a MATLAB toolbox for infinite integral 
% of products of two Bessel functions. To appear in ACM Transactions on 
% Mathematical Software. 


format long e
addpath('..','./BESSELINT')

%CC1 - 25
clear all
close all

rho_vec = [0.001 0.01 0.1 1 10 100 1000];
tau_vec = [0.0011 0.011 0.11 1.1 11 101 1001];

abserr = 1e-10;relerr = 1.0e-11;
u = 10+2i; a = 0; b = 0;type='JJ';
fx = @(x)x.*besselk(0,x.*u);
for i=1:length(rho_vec)
    for j=1:length(tau_vec)
        [results(i,j),rel_err(i,j), neval(i,j)]=IIPBF(fx,rho_vec(i),tau_vec(j),a,b,abserr,relerr,type);
        answer=1/sqrt((u^2 + rho_vec(i)^2+tau_vec(j)^2)^2-4*rho_vec(i)^2*tau_vec(i)^2);
        act_err(i,j) = (abs(results(i,j) - answer));
        act_err_1(i,j) = act_err(i,j)/(10.^floor(log10(act_err(i,j))));
        act_err_2(i,j) = floor(log10(act_err(i,j)));
        if (act_err(i,j)==0)
            act_err_1(i,j) = 0;
            act_err_2(i,j) = 0;
        end
    end
end

fprintf('\n\\multicolumn{9}{c}{Case 1A} \\\\\n');
fprintf('       &        &               \\multicolumn{7}{c}{$\\rho$} \\\\\n');
fprintf('       &        & %s & %s & %s & %s & %s & %s & %s \\\\ \\cline{3-9}\n', num2str(rho_vec(1)), num2str(rho_vec(2)), num2str(rho_vec(3)), num2str(rho_vec(4)), num2str(rho_vec(5)), num2str(rho_vec(6)), num2str(rho_vec(7)));
fprintf('       & %s & %1.2f %s & %1.2f %s & %1.2f %s & %1.2f %s & %1.2f %s & %1.2f %s & %1.2f %s \\\\\n', num2str(tau_vec(1)), act_err_1(1,1), ['x$10^{' num2str(act_err_2(1,1)) '}$' ],  act_err_1(2,1), ['x$10^{' num2str(act_err_2(2,1)) '}$' ], act_err_1(3,1), ['x$10^{' num2str(act_err_2(3,1)) '}$' ], act_err_1(4,1), ['x$10^{' num2str(act_err_2(4,1)) '}$' ], act_err_1(5,1), ['x$10^{' num2str(act_err_2(5,1)) '}$' ], act_err_1(6,1), ['x$10^{' num2str(act_err_2(6,1)) '}$' ], act_err_1(7,1), ['x$10^{' num2str(act_err_2(7,1)) '}$' ]);
fprintf('       & %s & %1.2f %s & %1.2f %s & %1.2f %s & %1.2f %s & %1.2f %s & %1.2f %s & %1.2f %s \\\\\n', num2str(tau_vec(2)), act_err_1(1,2), ['x$10^{' num2str(act_err_2(1,2)) '}$' ],  act_err_1(2,2), ['x$10^{' num2str(act_err_2(2,2)) '}$' ], act_err_1(3,2), ['x$10^{' num2str(act_err_2(3,2)) '}$' ], act_err_1(4,2), ['x$10^{' num2str(act_err_2(4,2)) '}$' ], act_err_1(5,2), ['x$10^{' num2str(act_err_2(5,2)) '}$' ], act_err_1(6,2), ['x$10^{' num2str(act_err_2(6,2)) '}$' ], act_err_1(7,2), ['x$10^{' num2str(act_err_2(7,2)) '}$' ]);
fprintf('       & %s & %1.2f %s & %1.2f %s & %1.2f %s & %1.2f %s & %1.2f %s & %1.2f %s & %1.2f %s \\\\\n', num2str(tau_vec(3)), act_err_1(1,3), ['x$10^{' num2str(act_err_2(1,3)) '}$' ],  act_err_1(2,3), ['x$10^{' num2str(act_err_2(2,3)) '}$' ], act_err_1(3,3), ['x$10^{' num2str(act_err_2(3,3)) '}$' ], act_err_1(4,3), ['x$10^{' num2str(act_err_2(4,3)) '}$' ], act_err_1(5,3), ['x$10^{' num2str(act_err_2(5,3)) '}$' ], act_err_1(6,3), ['x$10^{' num2str(act_err_2(6,3)) '}$' ], act_err_1(7,3), ['x$10^{' num2str(act_err_2(7,3)) '}$' ]);
fprintf('$\\tau$ & %s & %1.2f %s & %1.2f %s & %1.2f %s & %1.2f %s & %1.2f %s & %1.2f %s & %1.2f %s \\\\\n', num2str(tau_vec(4)), act_err_1(1,4), ['x$10^{' num2str(act_err_2(1,4)) '}$' ],  act_err_1(2,4), ['x$10^{' num2str(act_err_2(2,4)) '}$' ], act_err_1(3,4), ['x$10^{' num2str(act_err_2(3,4)) '}$' ], act_err_1(4,4), ['x$10^{' num2str(act_err_2(4,4)) '}$' ], act_err_1(5,4), ['x$10^{' num2str(act_err_2(5,4)) '}$' ], act_err_1(6,4), ['x$10^{' num2str(act_err_2(6,4)) '}$' ], act_err_1(7,4), ['x$10^{' num2str(act_err_2(7,4)) '}$' ]);
fprintf('       & %s & %1.2f %s & %1.2f %s & %1.2f %s & %1.2f %s & %1.2f %s & %1.2f %s & %1.2f %s \\\\\n', num2str(tau_vec(5)), act_err_1(1,5), ['x$10^{' num2str(act_err_2(1,5)) '}$' ],  act_err_1(2,5), ['x$10^{' num2str(act_err_2(2,5)) '}$' ], act_err_1(3,5), ['x$10^{' num2str(act_err_2(3,5)) '}$' ], act_err_1(4,5), ['x$10^{' num2str(act_err_2(4,5)) '}$' ], act_err_1(5,5), ['x$10^{' num2str(act_err_2(5,5)) '}$' ], act_err_1(6,5), ['x$10^{' num2str(act_err_2(6,5)) '}$' ], act_err_1(7,5), ['x$10^{' num2str(act_err_2(7,5)) '}$' ]);
fprintf('       & %s & %1.2f %s & %1.2f %s & %1.2f %s & %1.2f %s & %1.2f %s & %1.2f %s & %1.2f %s \\\\\n', num2str(tau_vec(6)), act_err_1(1,6), ['x$10^{' num2str(act_err_2(1,6)) '}$' ],  act_err_1(2,6), ['x$10^{' num2str(act_err_2(2,6)) '}$' ], act_err_1(3,6), ['x$10^{' num2str(act_err_2(3,6)) '}$' ], act_err_1(4,6), ['x$10^{' num2str(act_err_2(4,6)) '}$' ], act_err_1(5,6), ['x$10^{' num2str(act_err_2(5,6)) '}$' ], act_err_1(6,6), ['x$10^{' num2str(act_err_2(6,6)) '}$' ], act_err_1(7,6), ['x$10^{' num2str(act_err_2(7,6)) '}$' ]);
fprintf('       & %s & %1.2f %s & %1.2f %s & %1.2f %s & %1.2f %s & %1.2f %s & %1.2f %s & %1.2f %s \\\\ \\hline\n', num2str(tau_vec(7)), act_err_1(1,7), ['x$10^{' num2str(act_err_2(1,7)) '}$' ],  act_err_1(2,7), ['x$10^{' num2str(act_err_2(2,7)) '}$' ], act_err_1(3,7), ['x$10^{' num2str(act_err_2(3,7)) '}$' ], act_err_1(4,7), ['x$10^{' num2str(act_err_2(4,7)) '}$' ], act_err_1(5,7), ['x$10^{' num2str(act_err_2(5,7)) '}$' ], act_err_1(6,7), ['x$10^{' num2str(act_err_2(6,7)) '}$' ], act_err_1(7,7), ['x$10^{' num2str(act_err_2(7,7)) '}$' ]);

u = 2+2i; a = 0; b = 0;type='JJ';
fx = @(x)x.*besselk(0,x.*u);
for i=1:length(rho_vec)
    for j=1:length(tau_vec)
        [results(i,j),rel_err(i,j), neval(i,j)]=IIPBF(fx,rho_vec(i),tau_vec(j),a,b,abserr,relerr,type);
        answer=1/sqrt((u^2 + rho_vec(i)^2+tau_vec(j)^2)^2-4*rho_vec(i)^2*tau_vec(i)^2);
        act_err(i,j) = (abs(results(i,j) - answer));
        act_err_1(i,j) = act_err(i,j)/(10.^floor(log10(act_err(i,j))));
        act_err_2(i,j) = floor(log10(act_err(i,j)));
        if (act_err(i,j)==0)
            act_err_1(i,j) = 0;
            act_err_2(i,j) = 0;
        end
    end
end

fprintf('\n\\multicolumn{9}{c}{Case 1B} \\\\\n');
fprintf('       &        &               \\multicolumn{7}{c}{$\\rho$} \\\\\n');
fprintf('       &        & %s & %s & %s & %s & %s & %s & %s \\\\ \\cline{3-9}\n', num2str(rho_vec(1)), num2str(rho_vec(2)), num2str(rho_vec(3)), num2str(rho_vec(4)), num2str(rho_vec(5)), num2str(rho_vec(6)), num2str(rho_vec(7)));
fprintf('       & %s & %1.2f %s & %1.2f %s & %1.2f %s & %1.2f %s & %1.2f %s & %1.2f %s & %1.2f %s \\\\\n', num2str(tau_vec(1)), act_err_1(1,1), ['x$10^{' num2str(act_err_2(1,1)) '}$' ],  act_err_1(2,1), ['x$10^{' num2str(act_err_2(2,1)) '}$' ], act_err_1(3,1), ['x$10^{' num2str(act_err_2(3,1)) '}$' ], act_err_1(4,1), ['x$10^{' num2str(act_err_2(4,1)) '}$' ], act_err_1(5,1), ['x$10^{' num2str(act_err_2(5,1)) '}$' ], act_err_1(6,1), ['x$10^{' num2str(act_err_2(6,1)) '}$' ], act_err_1(7,1), ['x$10^{' num2str(act_err_2(7,1)) '}$' ]);
fprintf('       & %s & %1.2f %s & %1.2f %s & %1.2f %s & %1.2f %s & %1.2f %s & %1.2f %s & %1.2f %s \\\\\n', num2str(tau_vec(2)), act_err_1(1,2), ['x$10^{' num2str(act_err_2(1,2)) '}$' ],  act_err_1(2,2), ['x$10^{' num2str(act_err_2(2,2)) '}$' ], act_err_1(3,2), ['x$10^{' num2str(act_err_2(3,2)) '}$' ], act_err_1(4,2), ['x$10^{' num2str(act_err_2(4,2)) '}$' ], act_err_1(5,2), ['x$10^{' num2str(act_err_2(5,2)) '}$' ], act_err_1(6,2), ['x$10^{' num2str(act_err_2(6,2)) '}$' ], act_err_1(7,2), ['x$10^{' num2str(act_err_2(7,2)) '}$' ]);
fprintf('       & %s & %1.2f %s & %1.2f %s & %1.2f %s & %1.2f %s & %1.2f %s & %1.2f %s & %1.2f %s \\\\\n', num2str(tau_vec(3)), act_err_1(1,3), ['x$10^{' num2str(act_err_2(1,3)) '}$' ],  act_err_1(2,3), ['x$10^{' num2str(act_err_2(2,3)) '}$' ], act_err_1(3,3), ['x$10^{' num2str(act_err_2(3,3)) '}$' ], act_err_1(4,3), ['x$10^{' num2str(act_err_2(4,3)) '}$' ], act_err_1(5,3), ['x$10^{' num2str(act_err_2(5,3)) '}$' ], act_err_1(6,3), ['x$10^{' num2str(act_err_2(6,3)) '}$' ], act_err_1(7,3), ['x$10^{' num2str(act_err_2(7,3)) '}$' ]);
fprintf('$\\tau$ & %s & %1.2f %s & %1.2f %s & %1.2f %s & %1.2f %s & %1.2f %s & %1.2f %s & %1.2f %s \\\\\n', num2str(tau_vec(4)), act_err_1(1,4), ['x$10^{' num2str(act_err_2(1,4)) '}$' ],  act_err_1(2,4), ['x$10^{' num2str(act_err_2(2,4)) '}$' ], act_err_1(3,4), ['x$10^{' num2str(act_err_2(3,4)) '}$' ], act_err_1(4,4), ['x$10^{' num2str(act_err_2(4,4)) '}$' ], act_err_1(5,4), ['x$10^{' num2str(act_err_2(5,4)) '}$' ], act_err_1(6,4), ['x$10^{' num2str(act_err_2(6,4)) '}$' ], act_err_1(7,4), ['x$10^{' num2str(act_err_2(7,4)) '}$' ]);
fprintf('       & %s & %1.2f %s & %1.2f %s & %1.2f %s & %1.2f %s & %1.2f %s & %1.2f %s & %1.2f %s \\\\\n', num2str(tau_vec(5)), act_err_1(1,5), ['x$10^{' num2str(act_err_2(1,5)) '}$' ],  act_err_1(2,5), ['x$10^{' num2str(act_err_2(2,5)) '}$' ], act_err_1(3,5), ['x$10^{' num2str(act_err_2(3,5)) '}$' ], act_err_1(4,5), ['x$10^{' num2str(act_err_2(4,5)) '}$' ], act_err_1(5,5), ['x$10^{' num2str(act_err_2(5,5)) '}$' ], act_err_1(6,5), ['x$10^{' num2str(act_err_2(6,5)) '}$' ], act_err_1(7,5), ['x$10^{' num2str(act_err_2(7,5)) '}$' ]);
fprintf('       & %s & %1.2f %s & %1.2f %s & %1.2f %s & %1.2f %s & %1.2f %s & %1.2f %s & %1.2f %s \\\\\n', num2str(tau_vec(6)), act_err_1(1,6), ['x$10^{' num2str(act_err_2(1,6)) '}$' ],  act_err_1(2,6), ['x$10^{' num2str(act_err_2(2,6)) '}$' ], act_err_1(3,6), ['x$10^{' num2str(act_err_2(3,6)) '}$' ], act_err_1(4,6), ['x$10^{' num2str(act_err_2(4,6)) '}$' ], act_err_1(5,6), ['x$10^{' num2str(act_err_2(5,6)) '}$' ], act_err_1(6,6), ['x$10^{' num2str(act_err_2(6,6)) '}$' ], act_err_1(7,6), ['x$10^{' num2str(act_err_2(7,6)) '}$' ]);
fprintf('       & %s & %1.2f %s & %1.2f %s & %1.2f %s & %1.2f %s & %1.2f %s & %1.2f %s & %1.2f %s \\\\ \\hline\n', num2str(tau_vec(7)), act_err_1(1,7), ['x$10^{' num2str(act_err_2(1,7)) '}$' ],  act_err_1(2,7), ['x$10^{' num2str(act_err_2(2,7)) '}$' ], act_err_1(3,7), ['x$10^{' num2str(act_err_2(3,7)) '}$' ], act_err_1(4,7), ['x$10^{' num2str(act_err_2(4,7)) '}$' ], act_err_1(5,7), ['x$10^{' num2str(act_err_2(5,7)) '}$' ], act_err_1(6,7), ['x$10^{' num2str(act_err_2(6,7)) '}$' ], act_err_1(7,7), ['x$10^{' num2str(act_err_2(7,7)) '}$' ]);

u = 2+10i; a = 0; b = 0;type='JJ';
fx = @(x)x.*besselk(0,x.*u);
for i=1:length(rho_vec)
    for j=1:length(tau_vec)
        [results(i,j),rel_err(i,j), neval(i,j)]=IIPBF(fx,rho_vec(i),tau_vec(j),a,b,abserr,relerr,type);
        answer=1/sqrt((u^2 + rho_vec(i)^2+tau_vec(j)^2)^2-4*rho_vec(i)^2*tau_vec(i)^2);
        act_err(i,j) = (abs(results(i,j) - answer));
        act_err_1(i,j) = act_err(i,j)/(10.^floor(log10(act_err(i,j))));
        act_err_2(i,j) = floor(log10(act_err(i,j)));
        if (act_err(i,j)==0)
            act_err_1(i,j) = 0;
            act_err_2(i,j) = 0;
        end
    end
end

fprintf('\n\\multicolumn{9}{c}{Case 1C} \\\\\n');
fprintf('       &        &               \\multicolumn{7}{c}{$\\rho$} \\\\\n');
fprintf('       &        & %s & %s & %s & %s & %s & %s & %s \\\\ \\cline{3-9}\n', num2str(rho_vec(1)), num2str(rho_vec(2)), num2str(rho_vec(3)), num2str(rho_vec(4)), num2str(rho_vec(5)), num2str(rho_vec(6)), num2str(rho_vec(7)));
fprintf('       & %s & %1.2f %s & %1.2f %s & %1.2f %s & %1.2f %s & %1.2f %s & %1.2f %s & %1.2f %s \\\\\n', num2str(tau_vec(1)), act_err_1(1,1), ['x$10^{' num2str(act_err_2(1,1)) '}$' ],  act_err_1(2,1), ['x$10^{' num2str(act_err_2(2,1)) '}$' ], act_err_1(3,1), ['x$10^{' num2str(act_err_2(3,1)) '}$' ], act_err_1(4,1), ['x$10^{' num2str(act_err_2(4,1)) '}$' ], act_err_1(5,1), ['x$10^{' num2str(act_err_2(5,1)) '}$' ], act_err_1(6,1), ['x$10^{' num2str(act_err_2(6,1)) '}$' ], act_err_1(7,1), ['x$10^{' num2str(act_err_2(7,1)) '}$' ]);
fprintf('       & %s & %1.2f %s & %1.2f %s & %1.2f %s & %1.2f %s & %1.2f %s & %1.2f %s & %1.2f %s \\\\\n', num2str(tau_vec(2)), act_err_1(1,2), ['x$10^{' num2str(act_err_2(1,2)) '}$' ],  act_err_1(2,2), ['x$10^{' num2str(act_err_2(2,2)) '}$' ], act_err_1(3,2), ['x$10^{' num2str(act_err_2(3,2)) '}$' ], act_err_1(4,2), ['x$10^{' num2str(act_err_2(4,2)) '}$' ], act_err_1(5,2), ['x$10^{' num2str(act_err_2(5,2)) '}$' ], act_err_1(6,2), ['x$10^{' num2str(act_err_2(6,2)) '}$' ], act_err_1(7,2), ['x$10^{' num2str(act_err_2(7,2)) '}$' ]);
fprintf('       & %s & %1.2f %s & %1.2f %s & %1.2f %s & %1.2f %s & %1.2f %s & %1.2f %s & %1.2f %s \\\\\n', num2str(tau_vec(3)), act_err_1(1,3), ['x$10^{' num2str(act_err_2(1,3)) '}$' ],  act_err_1(2,3), ['x$10^{' num2str(act_err_2(2,3)) '}$' ], act_err_1(3,3), ['x$10^{' num2str(act_err_2(3,3)) '}$' ], act_err_1(4,3), ['x$10^{' num2str(act_err_2(4,3)) '}$' ], act_err_1(5,3), ['x$10^{' num2str(act_err_2(5,3)) '}$' ], act_err_1(6,3), ['x$10^{' num2str(act_err_2(6,3)) '}$' ], act_err_1(7,3), ['x$10^{' num2str(act_err_2(7,3)) '}$' ]);
fprintf('$\\tau$ & %s & %1.2f %s & %1.2f %s & %1.2f %s & %1.2f %s & %1.2f %s & %1.2f %s & %1.2f %s \\\\\n', num2str(tau_vec(4)), act_err_1(1,4), ['x$10^{' num2str(act_err_2(1,4)) '}$' ],  act_err_1(2,4), ['x$10^{' num2str(act_err_2(2,4)) '}$' ], act_err_1(3,4), ['x$10^{' num2str(act_err_2(3,4)) '}$' ], act_err_1(4,4), ['x$10^{' num2str(act_err_2(4,4)) '}$' ], act_err_1(5,4), ['x$10^{' num2str(act_err_2(5,4)) '}$' ], act_err_1(6,4), ['x$10^{' num2str(act_err_2(6,4)) '}$' ], act_err_1(7,4), ['x$10^{' num2str(act_err_2(7,4)) '}$' ]);
fprintf('       & %s & %1.2f %s & %1.2f %s & %1.2f %s & %1.2f %s & %1.2f %s & %1.2f %s & %1.2f %s \\\\\n', num2str(tau_vec(5)), act_err_1(1,5), ['x$10^{' num2str(act_err_2(1,5)) '}$' ],  act_err_1(2,5), ['x$10^{' num2str(act_err_2(2,5)) '}$' ], act_err_1(3,5), ['x$10^{' num2str(act_err_2(3,5)) '}$' ], act_err_1(4,5), ['x$10^{' num2str(act_err_2(4,5)) '}$' ], act_err_1(5,5), ['x$10^{' num2str(act_err_2(5,5)) '}$' ], act_err_1(6,5), ['x$10^{' num2str(act_err_2(6,5)) '}$' ], act_err_1(7,5), ['x$10^{' num2str(act_err_2(7,5)) '}$' ]);
fprintf('       & %s & %1.2f %s & %1.2f %s & %1.2f %s & %1.2f %s & %1.2f %s & %1.2f %s & %1.2f %s \\\\\n', num2str(tau_vec(6)), act_err_1(1,6), ['x$10^{' num2str(act_err_2(1,6)) '}$' ],  act_err_1(2,6), ['x$10^{' num2str(act_err_2(2,6)) '}$' ], act_err_1(3,6), ['x$10^{' num2str(act_err_2(3,6)) '}$' ], act_err_1(4,6), ['x$10^{' num2str(act_err_2(4,6)) '}$' ], act_err_1(5,6), ['x$10^{' num2str(act_err_2(5,6)) '}$' ], act_err_1(6,6), ['x$10^{' num2str(act_err_2(6,6)) '}$' ], act_err_1(7,6), ['x$10^{' num2str(act_err_2(7,6)) '}$' ]);
fprintf('       & %s & %1.2f %s & %1.2f %s & %1.2f %s & %1.2f %s & %1.2f %s & %1.2f %s & %1.2f %s \\\\ \\hline\n', num2str(tau_vec(7)), act_err_1(1,7), ['x$10^{' num2str(act_err_2(1,7)) '}$' ],  act_err_1(2,7), ['x$10^{' num2str(act_err_2(2,7)) '}$' ], act_err_1(3,7), ['x$10^{' num2str(act_err_2(3,7)) '}$' ], act_err_1(4,7), ['x$10^{' num2str(act_err_2(4,7)) '}$' ], act_err_1(5,7), ['x$10^{' num2str(act_err_2(5,7)) '}$' ], act_err_1(6,7), ['x$10^{' num2str(act_err_2(6,7)) '}$' ], act_err_1(7,7), ['x$10^{' num2str(act_err_2(7,7)) '}$' ]);

%CC2 - 26 (Gradshteyn & Ryzhik 6.633.3)
clear all; close all;
format long e;
type = 'JY';u = 10;b=1/3;rho = 1;tau = 1;
[res2,rel2,evals2]=IIPBF(@(x)power(x,2*b+1).*exp(-u.*x.*x),rho,tau,b,b,1.e-14,1.e-13,type);
f2 = -power(u,-(3*b+1)/2)*exp(-1/(2*u))*whittakerW(b/2,b/2,1/u)/(2*sqrt(pi));
diff2 = abs(res2-f2);
fprintf('2A & %10.6e & %10.6e \n', rel2, diff2);
clear all; close all;
format long e;
type = 'JY';u = 10;b=1/2;rho = 1;tau = 1;
[res2,rel2,evals2]=IIPBF(@(x)power(x,2*b+1).*exp(-u.*x.*x),rho,tau,b,b,1.e-14,1.e-13,type);
f2 = -power(u,-(3*b+1)/2)*exp(-1/(2*u))*whittakerW(b/2,b/2,1/u)/(2*sqrt(pi));
diff2 = abs(res2-f2);
fprintf('2B & %10.6e & %10.6e \n', rel2, diff2);

%CC3 -27 (Eq 31. McPhedran 1992) 
clear all; close all;
format long e;
type = 'JY';u = 10;b=1/3;rho = 1;tau = 1;
[res3,rel3,evals3]=IIPBF(@(x)x.*exp(-u.*x.*x),rho,tau,b,b,1.e-14,1.e-13,type);
LukeAssocBessel1 = -exp(1/(2*u))*hypergeom([1,1/2],[1-b,b+1],-1/u)/power(b,2);
f3 = (exp(-1/(2*u))/(2*pi*u))*(pi*cot(pi*b)*besseli(b,1/(2*u))+b*LukeAssocBessel1);
diff3 = abs(res3-f3);
fprintf('3A & %10.6e & %10.6e \n', rel3, diff3);
clear all; close all;
format long e;
type = 'JY';u = 10;b=1/2;rho = 1;tau = 1;
[res3,rel3,evals3]=IIPBF(@(x)x.*exp(-u.*x.*x),rho,tau,b,b,1.e-14,1.e-13,type);
LukeAssocBessel1 = -exp(1/(2*u))*hypergeom([1,1/2],[1-b,b+1],-1/u)/power(b,2);
f3 = (exp(-1/(2*u))/(2*pi*u))*(pi*cot(pi*b)*besseli(b,1/(2*u))+b*LukeAssocBessel1);
diff3 = abs(res3-f3);
fprintf('3B & %10.6e & %10.6e \n', rel3, diff3);
clear all; close all;
format long e;
type = 'JY';u = 10;b=3/2;rho = 1;tau = 1;
[res3,rel3,evals3]=IIPBF(@(x)x.*exp(-u.*x.*x),rho,tau,b,b,1.e-14,1.e-13,type);
LukeAssocBessel1 = -exp(1/(2*u))*hypergeom([1,1/2],[1-b,b+1],-1/u)/power(b,2);
f3 = (exp(-1/(2*u))/(2*pi*u))*(pi*cot(pi*b)*besseli(b,1/(2*u))+b*LukeAssocBessel1);
diff3 = abs(res3-f3);
fprintf('3C & %10.6e & %10.6e \n', rel3, diff3);

%CC4 - 29 (Eq. 38 - McPhedran 1992)
clear all; close all;
format long e;
type = 'JJ';u = 10;b=1/3;rho = 1;tau = 1;
[res4,rel4,evals4]=IIPBF(@(x)x.*exp(-u.*x.*x),rho,tau,b,-b,1.e-14,1.e-13,type);
LukeAssocBessel1 = -exp(1/(2*u))*hypergeom([1,1/2],[1-b,b+1],-1/u)/power(b,2);
f4 = -(exp(-1/(2*u))/(2*pi*u))*(b*sin(pi*b)*LukeAssocBessel1);
diff4 = abs(res4-f4);
fprintf('4A & %10.6e & %10.6e \n', rel4, diff4);
clear all; close all;
format long e;
type = 'JJ';u = 10;b=1/2;rho = 1;tau = 1;
[res4,rel4,evals4]=IIPBF(@(x)x.*exp(-u.*x.*x),rho,tau,b,-b,1.e-14,1.e-13,type);
LukeAssocBessel1 = -exp(1/(2*u))*hypergeom([1,1/2],[1-b,b+1],-1/u)/power(b,2);
f4 = -(exp(-1/(2*u))/(2*pi*u))*(b*sin(pi*b)*LukeAssocBessel1);
diff4 = abs(res4-f4);
fprintf('4B & %10.6e & %10.6e \n', rel4, diff4);
clear all; close all;
format long e;
type = 'JJ';u = 10;b=3/2;rho = 1;tau = 1;
[res4,rel4,evals4]=IIPBF(@(x)x.*exp(-u.*x.*x),rho,tau,b,-b,1.e-14,1.e-13,type);
LukeAssocBessel1 = -exp(1/(2*u))*hypergeom([1,1/2],[1-b,b+1],-1/u)/power(b,2);
f4 = -(exp(-1/(2*u))/(2*pi*u))*(b*sin(pi*b)*LukeAssocBessel1);
diff4 = abs(res4-f4);
fprintf('4C & %10.6e & %10.6e \n', rel4, diff4);

%CC5 - 30 (Eq. 39 - McPhedran 1992)
clear all; close all;
format long e;
type = 'JY';u = 10;b=1/3;rho = 1;tau = 1;
[res5,rel5,evals5]=IIPBF(@(x)x.*exp(-u.*x.*x),rho,tau,b,-b,1.e-14,1.e-13,type);
LukeAssocBessel1 = -exp(1/(2*u))*hypergeom([1,1/2],[1-b,b+1],-1/u)/power(b,2);
f5 = (exp(-1/(2*u))/(2*pi*u))*(pi*csc(pi*b)*besseli(b,1/(2*u))+b*cos(pi*b)*LukeAssocBessel1);
diff5 = abs(res5-f5);
fprintf('5A & %10.6e & %10.6e \n', rel5, diff5);
clear all; close all;
format long e;
type = 'JY';u = 10;b=1/2;rho = 1;tau = 1;
[res5,rel5,evals5]=IIPBF(@(x)x.*exp(-u.*x.*x),rho,tau,b,-b,1.e-14,1.e-13,type);
LukeAssocBessel1 = -exp(1/(2*u))*hypergeom([1,1/2],[1-b,b+1],-1/u)/power(b,2);
f5 = (exp(-1/(2*u))/(2*pi*u))*(pi*csc(pi*b)*besseli(b,1/(2*u))+b*cos(pi*b)*LukeAssocBessel1);
diff5 = abs(res5-f5);
fprintf('5B & %10.6e & %10.6e \n', rel5, diff5);
clear all; close all;
format long e;
type = 'JY';u = 10;b=3/2;rho = 1;tau = 1;
[res5,rel5,evals5]=IIPBF(@(x)x.*exp(-u.*x.*x),rho,tau,b,-b,1.e-14,1.e-13,type);
LukeAssocBessel1 = -exp(1/(2*u))*hypergeom([1,1/2],[1-b,b+1],-1/u)/power(b,2);
f5 = (exp(-1/(2*u))/(2*pi*u))*(pi*csc(pi*b)*besseli(b,1/(2*u))+b*cos(pi*b)*LukeAssocBessel1);
diff5 = abs(res5-f5);
fprintf('5C & %10.6e & %10.6e \n', rel5, diff5);

%CC6 - 31 (Eq. Eq. 40 - McPhedran 1992)
clear all; close all;
format long e;
type = 'YY';u = 10;b=1/3;rho = 1;tau = 1;
[res6,rel6,evals6]=IIPBF(@(x)x.*exp(-u.*x.*x),rho,tau,b,-b,1.e-14,1.e-13,type);
LukeAssocBessel1 = -exp(1/(2*u))*hypergeom([1,1/2],[1-b,b+1],-1/u)/power(b,2);
f6 = (exp(-1/(2*u))/(2*pi*u))*(pi*cos(pi*b)*(besseli(b,1/(2*u))+besseli(-b,1/(2*u)))+(1+power(cos(b*pi),2))*b*sin(pi*b)*LukeAssocBessel1)/power(sin(b*pi),2);
diff6 = abs(res6-f6);
fprintf('6A & %10.6e & %10.6e \n', rel6, diff6);
clear all; close all;
format long e;
type = 'YY';u = 10;b=1/2;rho = 1;tau = 1;
[res6,rel6,evals6]=IIPBF(@(x)x.*exp(-u.*x.*x),rho,tau,b,-b,1.e-14,1.e-13,type);
LukeAssocBessel1 = -exp(1/(2*u))*hypergeom([1,1/2],[1-b,b+1],-1/u)/power(b,2);
f6 = (exp(-1/(2*u))/(2*pi*u))*(pi*cos(pi*b)*(besseli(b,1/(2*u))+besseli(-b,1/(2*u)))+(1+power(cos(b*pi),2))*b*sin(pi*b)*LukeAssocBessel1)/power(sin(b*pi),2);
diff6 = abs(res6-f6);
fprintf('6B & %10.6e & %10.6e \n', rel6, diff6);

%CC7 - 32 (Eq. 41 - McPhedran 1992)
clear all; close all;
format long e;
type = 'YY';u = 10;b=1/3;rho = 1;tau = 1;
[res7,rel7,evals7]=IIPBF(@(x)x.*exp(-u.*x.*x),rho,tau,b,b,1.e-14,1.e-13,type);
LukeAssocBessel1 = -exp(1/(2*u))*hypergeom([1,1/2],[1-b,b+1],-1/u)/power(b,2);
f7 = (exp(-1/(2*u))/(2*pi*u))*(pi*power(cot(pi*b),2)*besseli(b,1/(2*u)) + pi*(1+power(cot(b*pi),2))*besseli(-b,1/(2*u))+ 2*b*cot(pi*b)*LukeAssocBessel1);
diff7 = abs(res7-f7);
fprintf('7A & %10.6e & %10.6e \n', rel7, diff7);
clear all; close all;
format long e;
type = 'YY';u = 10;b=1/2;rho = 1;tau = 1;
[res7,rel7,evals7]=IIPBF(@(x)x.*exp(-u.*x.*x),rho,tau,b,b,1.e-14,1.e-13,type);
LukeAssocBessel1 = -exp(1/(2*u))*hypergeom([1,1/2],[1-b,b+1],-1/u)/power(b,2);
f7 = (exp(-1/(2*u))/(2*pi*u))*(pi*power(cot(pi*b),2)*besseli(b,1/(2*u)) + pi*(1+power(cot(b*pi),2))*besseli(-b,1/(2*u))+ 2*b*cot(pi*b)*LukeAssocBessel1);
diff7 = abs(res7-f7);
fprintf('7B & %10.6e & %10.6e \n', rel7, diff7);

%CC8 33 (Eq. 7 McPhedran 1992)
clear all; close all;
format long e;
type = 'JY';u = 10;b=0;rho = 1;tau = 1;
[res8,rel8,evals8]=IIPBF(@(x)x.*exp(-u.*x.*x),rho,tau,b,b,1.e-14,1.e-13,type);
f8 = -exp(-1/(2*u))*besselk(0,1/(2*u))/(2*pi*u);
diff8 = abs(res8-f8);
fprintf('8 & %10.6e & %10.6e \n', rel8, diff8);

%CC9 34 (Eq. 34 McPhedran 1992)
clear all; close all;
format long e;
type = 'JY';u = 0.1;b=10;rho = 1;tau = 1;
[res9,rel9,evals9]=IIPBF(@(x)x.*exp(-u.*x.*x),rho,tau,b,b,1.e-14,1.e-13,type);
z = -1/(2*u);
LukeAssocBessel2 = (exp(-z)/z)*hypergeom([1,b+1,1-b],3/2,-1/(2*z));
f9 = z*exp(z)*(power(-1,b)*besselk(b,-z)-b*LukeAssocBessel2)/pi;
diff9 = abs(res9-f9);
fprintf('9 & %10.6e & %10.6e \n', rel9, diff9);
