function y = icdf_PSM( u, param, maxmin, du )
    arg = maxmin(1):du:maxmin(2);
    y1 = cdf_PSM( arg, param );
    y=interp1(y1,arg,u,'PCHIP','extrap');
return;
