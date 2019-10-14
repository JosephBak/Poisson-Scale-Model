function y = cdf_stdPSM_old(x, stdparam, M)
    if nargin == 2
        M = 10;
    end
    
    param = convertStdparam2Nonstdparam(stdparam);
    if length(stdparam)>=3
        y=cdf_PSM_old(x, [param, stdparam(3)], M);
    else
        y=cdf_PSM_old(x, param, M);
    end
 end