datastyle = 1;
method = 1;
    
    SymbolFactors = { '^GSPC' };
    yDataIndex = getYahooDailyData(SymbolFactors, '01/01/2013', '12/31/2015', 'mm/dd/yyyy');
    datalength = length(yDataIndex.(genvarname(SymbolFactors{1})).AdjClose);
    numoffactor = length(SymbolFactors);
    retIndex = ...
        log(yDataIndex.(genvarname(SymbolFactors{1})).AdjClose(2:datalength)./yDataIndex.(genvarname(SymbolFactors{1})).AdjClose(1:datalength-1));

    
    if datastyle == 1%std case
        stdret = (retIndex-mean(retIndex))/std(retIndex);
        inputdata = stdret;
    else %non-std case
        inputdata = retIndex;
    end
    
    [empCumDis,xi_empCumDis] = ksdensity(inputdata(:),'function','cdf');
    [empFreDis,xi] = ksdensity(inputdata(:));

    if datastyle == 1%std case
        if method == 1  %Curvefit cdf of stdPSM
            initial = [0.0, 1];
            [xlsq, resnorm] = lsqcurvefit( @(x,y)cdf_stdPSM(y, x), initial, xi_empCumDis, empCumDis, [ -0.5, 0.250001], [0.5, Inf]);
        elseif method == 2  %Curvefit pdf of stdPSM
            initial = [0.0, 1];
            [xlsq, resnorm] = lsqcurvefit( @(x,y)pdf_stdPSM(y, x), initial, xi, empFreDis, [ -0.5, 0.250001], [0.5, Inf]);
        elseif method == 3 % curve fit chf of stdPSM
            wk = 0.1:0.1:1;
            wk = wk(:);
            samplelogsymcf = log(samplechf(wk, stdret));
            initial = [0.0, 1];
            opt = optimoptions(@lsqnonlin, 'TolX', 1e-9, 'TolFun', 1e-9, 'MaXFunEvals', 1000);
            [xlsq, resnorm] = lsqnonlin(@(x)errstcf(x, samplelogsymcf, wk, @chf_stdPSM), initial, [-0.5,0.250001], [0.5,inf], opt);
        end
    else %non-std case
        initial = [0.0, 1, 1.4, 1];
        if method == 1  %Curvefit cdf of PSM
            [xlsq, resnorm] = lsqcurvefit( @(x,y)cdf_PSM(y, x), initial, xi_empCumDis, empCumDis, [-Inf, -inf, 0, 0], [Inf, Inf,Inf, Inf]);
        elseif method == 2  %Curvefit pdf of PSM
            [xlsq, resnorm] = lsqcurvefit( @(x,y)pdf_PSM(y, x), initial, xi, empFreDis, [-Inf, -inf, 0, 0], [Inf, Inf,Inf, Inf]);
        end
    end
    
    pdfPSM = pdf_stdPSM( xi, xlsq );  
    cdfPSM = cdf_stdPSM( xi, xlsq );  
    pdfNormal = pdf('norm', xi, 0, 1);
    cdfG = cdf('norm', xi_empCumDis, 0, 1);
    
    figure(1),plot(xi, empFreDis, 'k-', xi, pdfPSM, 'k--', xi, pdfNormal, 'k-.')
    figure(2),plot(xi, log(empFreDis), 'k', xi, log(pdfPSM), 'r', xi, log(pdfNormal), 'b') 
    figure(3),plot(xi_empCumDis, empCumDis, 'k', xi_empCumDis, cdfPSM, 'r', xi_empCumDis, cdfG, 'b') 

    [h,pPSM,ksstatPSM,cv] = kstest( inputdata(:),[xi_empCumDis(:),cdfPSM(:)] );
    [adPSM, ad2PSM] = ad2test_FFT( inputdata(:),[xi_empCumDis(:),cdfPSM(:)] );
    ad2PSM=ad2PSM*length(inputdata);
    pAD2PSM=1-AD(length(inputdata),ad2PSM);
    
    [h,pG,ksstatG,cv] = kstest(inputdata(:),[xi_empCumDis(:),cdfG(:)]);
    [adG, ad2G] = ad2test_FFT( inputdata(:),[xi_empCumDis(:),cdfG(:)] );
    ad2G=ad2G*length(inputdata);
    pAD2G=1-AD(length(inputdata),ad2G);
    
    figure(4)
    plot(-4:0.1:4, pdf_stdPSM( -4:0.1:4, xlsq )./pdf('norm', -4:0.1:4, 0, 1));
    figure(5)
    plot(-4:0.1:4, cdf_stdPSM( -4:0.1:4, xlsq )./cdf('norm', -4:0.1:4, 0, 1));

    
    return
