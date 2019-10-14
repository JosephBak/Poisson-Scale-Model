function price = CallPriceVectorPSM( param, xdata, S0, r, d )
    
    K = xdata(:,1);
    tau = xdata(:,2);
    
    price = zeros(size(K));
    difftau = diff(tau);
    selectedtau = [tau(difftau~=0); tau(length(tau))];
    indextau = [0; find(difftau~=0);length(tau)];
    price = zeros(size(K));
    for j = 1:length(selectedtau)
        inpK=K(indextau(j)+1:indextau(j+1));
        c = CallPutPriceCTS_FFT( inpK(:).', S0(1), r(1), d(1), selectedtau(j), param, 'C' ) ;
        price(indextau(j)+1:indextau(j+1)) = c;
    end
return