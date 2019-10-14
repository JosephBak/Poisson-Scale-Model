SymbolFactors = { '^GSPC' };
yDataIndex = getYahooDailyData(SymbolFactors, '01/01/2013', '12/31/2015', 'mm/dd/yyyy');
datalength = length(yDataIndex.(genvarname(SymbolFactors{1})).AdjClose);
numoffactor = length(SymbolFactors);
retIndex = ...
    log(yDataIndex.(genvarname(SymbolFactors{1})).AdjClose(2:datalength)./yDataIndex.(genvarname(SymbolFactors{1})).AdjClose(1:datalength-1));

    stdret = (retIndex-mean(retIndex))/std(retIndex);

    initial = [1, 1.4, 1, 0.0];

    [empCumDis,xi_empCumDis] = ksdensity(stdret(:),'function','cdf');
    [empFreDis,xi] = ksdensity(stdret(:));
    [xlsq, resnorm] = lsqcurvefit( @(x,y)cdf_PSM(y, x), initial, xi_empCumDis, empCumDis, [0, -Inf, 0, -inf], [Inf, Inf,Inf, Inf]);
    
    pdfPSM = pdf_PSM( xi, xlsq );  
    cdfPSM = cdf_PSM( xi, xlsq );  

    figure(1),plot(xi, empFreDis, 'k', xi, pdfPSM, 'b') 
    figure(2),plot(xi_empCumDis, empCumDis, 'k', xi_empCumDis, cdfPSM, 'b') 

    [h,pPSM,ksstatPSM,cv] = kstest(stdret(:),[xi_empCumDis(:),cdfPSM(:)]);
    [adPSM, ad2PSM] = ad2test_FFT( stdret(:),[xi_empCumDis(:),cdfPSM(:)] );
    ad2PSM=ad2PSM*length(stdret);
    pAD2PSM=1-AD(length(stdret),ad2PSM);
return
