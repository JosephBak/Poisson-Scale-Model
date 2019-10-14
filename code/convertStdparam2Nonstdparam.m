function convparam = convertStdparam2Nonstdparam(stdparam)
    gam   = stdparam(1);
    lamb   = stdparam(2);
    %when k=0
    sig  = sqrt((1-lamb*gam^2)/(lamb));
    mu = -gam*lamb;
    
    %when k =2
    %sig  = sqrt((1-lamb*gam^2)/(2+lamb));
    %mu = -2*gam-gam*lamb;
    convparam = [mu, gam, lamb, sig];
return;