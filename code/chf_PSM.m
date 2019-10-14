function phi = chf_PSM(u, param)

    mu   = param(1);
    gam  = param(2);
    lamb = param(3);
    sig  = param(4);

    if length(param)>=5
        dt = param(5);
    else
        dt = 1;
    end
    
    phi = exp(dt*(1i*u*mu+2*u*gam*1i-sig^2*u.^2-lamb*(1-exp(1i*u*gam-sig^2*u.^2/2))));

return;