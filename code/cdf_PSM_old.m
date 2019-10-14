function y = cdf_PSM_old(x, param, M)
    mu   = param(1);
    gam     = param(2);
    lamb = param(3);
    sig   = param(4);
    
    n = length(x);
    if nargin == 2
        M = 20;
    end
    J = 0:M;
    lambda_jfac = (lamb.^(J)./factorial(J)).*exp(-lamb);
    
    gamj = gam*(J+2);
    s = sig^2*(J+2);
    m = mu+gamj;
    
    y = zeros(length(x), length(J));
    for k = 1:length(x)
        y(k, :) = cdf('norm', x(k), m(:), sqrt(s(:))); 
    end
    
    y = y*lambda_jfac(:);
    [a,b] = size(x);
    if a<b
        y = y.';
    end
 end