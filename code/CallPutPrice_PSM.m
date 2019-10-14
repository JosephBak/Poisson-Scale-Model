function [ OptionPrice ] = CallPutPrice_PSM( K, S0, r, tau, param ,flag)
    
    %param(1) is the sigma.
    %param(2) is the lamda.
    syms k
    if (flag == 'C' || flag == 'c')
        OptionPrice = vpa(exp(-param(2)*tau)*symsum((S0*normcdf(d_1(K, S0, k, param)) ...
        - K*exp(-r*tau)*normcdf(d_2(K, S0, k, param)))*(((param(2)*tau)^k)/factorial(k)), k, 1, Inf), 10);
    %the precision is 10.
    elseif (flag == 'P' || flag == 'p')
        OptionPrice = vpa(exp(-param(2)*tau)*symsum((K*exp(-r*tau)*normcdf(-d_2(K, S0, k, param)) ...
            - S0*normcdf(-d_1(K, S0, k, param)))*(((param(2)*tau)^k)/factorial(k)),k, 1, Inf),10);
     %the precision is 10.
    else
        disp( 'Flag must be C or P' );
    end
     
return
function d1 = d_1(K, S0, n, param)
    d1 = (log(S0/K) + ((param(1)^2)/2)*n)/(param(1)*sqrt(n));
return
function d2 = d_2(K, S0, n, param)
    d2 = d_1(K, S0, n, param) - param(1)*sqrt(n);
return