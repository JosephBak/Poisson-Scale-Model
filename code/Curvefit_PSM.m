    SymbolFactors = { '^GSPC' };
    yDataIndex = getYahooDailyData(SymbolFactors, '01/01/2013', '12/31/2015', 'mm/dd/yyyy');
    datalength = length(yDataIndex.(genvarname(SymbolFactors{1})).AdjClose);
    numoffactor = length(SymbolFactors);
    retIndex = ...
        log(yDataIndex.(genvarname(SymbolFactors{1})).AdjClose(2:datalength)./yDataIndex.(genvarname(SymbolFactors{1})).AdjClose(1:datalength-1));

    stdret = (retIndex-mean(retIndex))/std(retIndex);
    inputdata = stdret;

    [empStdCumDis,xi_empStdCumDis] = ksdensity(stdret(:),'function','cdf');
    [empStdFreDis,xi_empStdFreDis] = ksdensity(stdret(:));
    initial = [0.0, 1];
    [xlsq, resnorm] = lsqcurvefit( @(x,y)cdf_stdPSM(y, x), initial, xi_empStdCumDis, empStdCumDis, [ -0.5, 0.250001], [0.5, Inf]);
       
    [empCumDis,xi_empCumDis] = ksdensity(retIndex(:),'function','cdf');
    [empFreDis,xi_empFreDis] = ksdensity(retIndex(:));

    conv_param = convertParamLinCom(xlsq, std(retIndex), mean(retIndex) )
    
    pdfPSM = pdf_PSM( xi_empFreDis, conv_param );  
    cdfPSM = cdf_PSM( xi_empCumDis, conv_param );  
    pdfNormal = pdf('norm', xi_empFreDis, mean(retIndex), std(retIndex));
    cdfG = cdf('norm', xi_empCumDis, mean(retIndex), std(retIndex));
    
    figure(1),plot(xi_empFreDis, empFreDis, 'k-', xi_empFreDis, pdfPSM, 'k--', xi_empFreDis, pdfNormal, 'k-.'), ...
        legend('Empirical distribution', 'PSM distribution', 'Gaussian distribution')
    figure(2),plot(xi_empFreDis, log(empFreDis), 'k', xi_empFreDis, log(pdfPSM), 'k--', xi_empFreDis, log(pdfNormal), 'k-.'), ...
        legend('Empirical distribution', 'PSM distribution', 'Gaussian distribution') 
    figure(3),plot(xi_empCumDis, empCumDis, 'k', xi_empCumDis, cdfPSM, 'k--', xi_empCumDis, cdfG, 'k-.') , ...
        legend('Empirical distribution', 'PSM distribution', 'Gaussian distribution')
    
    u = 0.0001:0.0001:0.9999;
    xmin = min(retIndex);
    xmax = max(retIndex);
    %QQPLOT PSM
    figure(4)
    ic = icdf_PSM( u, conv_param, [xmin, xmax], (xmax-xmin)/100 );
    h = qqplotchr(retIndex, ic);
    axis([xmin, xmax, xmin, xmax]), title('PSM');
    %QQPLOT PSM
    figure(5)
    icG = icdf( 'norm', u, mean(retIndex), std(retIndex) );
    h = qqplotchr(retIndex, icG);
    axis([xmin, xmax, xmin, xmax]), title('Gaussian');
    
    [h,pPSM,ksstatPSM,cv] = kstest( retIndex(:),[xi_empCumDis(:),cdfPSM(:)] );
    [adPSM, ad2PSM] = ad2test_FFT( retIndex(:),[xi_empCumDis(:),cdfPSM(:)] );
    ad2PSM=ad2PSM*length(retIndex);
    pAD2PSM=1-AD(length(retIndex),ad2PSM);
    
    cdfG = cdf('norm', xi_empStdCumDis, 0, 1);
    [h,pG,ksstatG,cv] = kstest(retIndex(:),[xi_empCumDis(:),cdfG(:)]);
    [adG, ad2G] = ad2test_FFT( retIndex(:),[xi_empCumDis(:),cdfG(:)] );
    ad2G=ad2G*length(retIndex);
    pAD2G=1-AD(length(retIndex),ad2G);
    
    figure(6)
    plot(xi_empFreDis, pdf_PSM( xi_empFreDis, conv_param )./pdf('norm', xi_empFreDis, mean(retIndex), std(retIndex)));
    figure(7)
    plot(xi_empCumDis, cdf_PSM( xi_empCumDis, conv_param )./cdf('norm', xi_empCumDis, mean(retIndex), std(retIndex)));
    
    
    