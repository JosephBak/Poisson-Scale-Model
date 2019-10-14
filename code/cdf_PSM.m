function y = cdf_PSM(x, param, M)
    mu   = param(1);
    gam     = param(2);
    lamb = param(3);
    sig   = param(4);
    
    n = length(x);
    if nargin == 2
        M = 20;
    end
    
    J = 0:M;
    gamj = gam*(repmat(J, n, 1)+2);
    sig2j = sig^2*(repmat(J, n, 1)+2);
        
    lambda_jfac = (repmat(lamb.^(J), n, 1)./repmat(factorial(J), n,1)).*exp(-lamb);
    X = repmat(x(:), 1,length(J));
    
    normalizesX = (X-mu-gamj)./sqrt(sig2j);
    y = cdf('norm', normalizesX, 0, 1); 
    
    y = y.*lambda_jfac;
    y = sum(y, 2);
    [a,b] = size(x);
    if a<b
        y = y.';
    end
 end