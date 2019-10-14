function y = cdf_stdPSM(x, stdparam, M)
    if nargin == 2
        M = 10;
    end
    
    param = convertStdparam2Nonstdparam(stdparam);
    if length(stdparam)>=3
        y=cdf_PSM(x, [param, stdparam(3)], M);
    else
        y=cdf_PSM(x, param, M);
    end
 end