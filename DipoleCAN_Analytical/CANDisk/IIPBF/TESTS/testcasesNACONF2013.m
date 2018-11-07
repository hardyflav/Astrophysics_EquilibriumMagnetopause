% Test cases 6 for complex-valued c, 33-35 for IIPBF 
% from Gradshteyn & Ryzhik 7th Volume
% called by the testcases_errorplotsNACONF2013.m script
% Geoffrey Gunter, 6/3/2013

% Reference:
% Ratnanather, J. T., Kim, J. H., Zhang, S., Davis, A. M. J., and Lucas, 
% S. K. 2013. Algorithm XXX: IIPBF, a MATLAB toolbox for infinite integral 
% of products of two Bessel functions. To appear in ACM Transactions on 
% Mathematical Software.

function [solutions,answers,est_errs,nevals,required] = testcasesNACONF2013(n,varargin)
addpath('..')
tic;format long e;warning off

switch n
    case 6
        %% Case 6
        %  G&R 6.522.5
        
        % f(x) = x*K_0(x*c)
        % B = J_0(rho*x)*J_0(tau*x)
        % [Re(c) > |Im(rho)|, tau > 0]
        c = varargin{1};
        rho = varargin{2};
        tau = varargin{3};
        
        fx = @(x) x.*besselk(0,c.*x);
        a = 0;
        b = 0;
        type = 'JJ';
        
        if real(c)>abs(imag(rho)) && tau>0
            answer = power(power(c^2+rho^2+tau^2,2) - 4*rho^2*tau^2,-(1/2));
        else
            error('Case 6: Re(c) > |Im(rho)|, tau > 0')
        end
        
    case 33
        %% Case 33
        %  G&R 6.541.1
        
        % f(x) = x/(x^2+c^2)
        % B = J_a(rho*x)*J_a(tau*x)
        % [Re(c) > 0, Re(a) > -1]
        c = varargin{1};
        a = varargin{2}; b = a;
        rho = varargin{3};
        tau = varargin{4};
        
        fx = @(x) x./(x.^2+c.^2);
        type = 'JJ';
        
        if real(c)>0 && real(a)>-1
            if tau<=rho
                answer = besseli(a,tau*c)*besselk(a,rho*c);
            elseif rho<tau
                answer = besseli(a,rho*c)*besselk(a,tau*c);
            end
        else
            error('Case 33: Re(c) > 0, Re(a) > -1')
        end
        
    case 34
        %% Case 34
        %  G&R 6.574.3
        
        % f(x) = x^-c
        % B = J_a(rho*x)*J_b(tau*x)
        % [Re(a+b-c+1) > 0, Re(c) > -1, tau < rho]
        c = varargin{1};
        a = varargin{2};
        b = varargin{3};
        rho = varargin{4};
        tau = varargin{5};
        
        fx = @(x) x.^(-c);
        type = 'JJ';
        
        if real(a+b-c+1)>0 && real(c)>-1 && rho>tau
            answer = (((tau.^b).*gamma((b+a-c+1)./2)) ...
                ./ ((2.^c).*(rho.^(b-c+1)).*gamma((a-b+c+1)./2).*gamma(b+1))) ...
                .* hypergeom([(a+b-c+1)./2, (-a+b-c+1)./2], b+1, ((tau.^2)./(rho.^2)));
        else
            error('Case 34: Re(a+b-c+1) > 0, Re(c) > -1, tau < rho')
        end
        
        
    case 35
        %% Case 35
        %  G&R 6.576.6
        
        % f(x) = x^-c
        % B = J_a(rho*x)*Y_b(tau*x)
        % [tau > rho, Re(c) > -1, Re(a-c+1+b) > 0, Re(a-c+1-b) > 0]
        c = varargin{1};
        a = varargin{2};
        b = varargin{3};
        rho = varargin{4};
        tau = varargin{5};

        fx = @(x) x.^(-c);
        type = 'JY';
        
        if tau>rho && real(c)>-1 && real(a-c+1+b)>0 && real(a-c+1-b)>0
            answer = (2./pi) .* sin(pi.*(a-b-c)./2) ...
                .* ((rho.^a.*gamma(.5-.5.*c+.5.*b+.5.*a).*gamma(.5-.5.*c-.5.*b+.5.*a)) ...
                ./ (2.^(c+1).*tau.^(-c+a+1).*gamma(a+1))) ...
                .* hypergeom([.5-.5.*c+.5.*b+.5.*a, .5-.5.*c-.5.*b+.5.*a], a+1, ((rho.^2)./(tau.^2)));
        else
            error('Case 35: tau > rho, Re(c) > -1, Re(a-c+1+b) > 0, Re(a-c+1-b) > 0')
        end

    otherwise
        error('Cases 6, 33-35 only')
end;

I=11;
if n==6
    if imag(c) > 3
        I=10;
    end
end

solutions = zeros(I,1);
answers = zeros(I,1);
est_errs = zeros(I,1);
nevals = zeros(I,1);
required = zeros(I,1);

for i = 1:I
    
    abserr = (1.0e-3)/power(10,i);
    relerr = (1.0e-11);
    
    [result,err,neval]=IIPBF(fx,rho,tau,a,b,abserr,relerr,type);
    solutions(i) = result;
    answers(i) = answer;
    est_errs(i) = err;
    nevals(i) = neval;
    required(i) = abserr;
    
    % fprintf('\nComputed ='), disp(result)
    % fprintf('Analytical ='), disp(answer)
    % fprintf('Absolute Error ='), disp(abs(result-answer))
    % fprintf('Estimated Error ='), disp(err)
    
end;

warning on;

toc;
