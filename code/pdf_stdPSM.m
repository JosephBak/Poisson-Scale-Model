function y = pdf_stdPSM(x, stdparam, M)
    if nargin == 2
        M = 10;
    end
    
    param = convertStdparam2Nonstdparam(stdparam)
    
    if length(stdparam)>=3
        y=pdf_PSM(x, [param, stdparam(3)], M);
    else
        y=pdf_PSM(x, param, M);
    end
 end