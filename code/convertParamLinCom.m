function conv_param = convertParamLinCom(stdparam, a, b )
    param = convertStdparam2Nonstdparam(stdparam);
    conv_param = [param(1)*a+b, a*param(2), param(3), abs(a)*param(4)];
end
