function y=chf_stdPSM(u, stdparam)
    param = convertStdparam2Nonstdparam(stdparam);
    if length(stdparam)>=3
        y=chf_PSM(u, [param, stdparam(3)]);
    else
        y=chf_PSM(u, param);
    end
return;
    