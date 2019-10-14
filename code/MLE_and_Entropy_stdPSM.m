SymbolFactors = { '^GSPC' };
yDataIndex = getYahooDailyData(SymbolFactors, '01/01/2013', '12/31/2015', 'mm/dd/yyyy');
datalength = length(yDataIndex.(genvarname(SymbolFactors{1})).AdjClose);
numoffactor = length(SymbolFactors);
retIndex = ...
    log(yDataIndex.(genvarname(SymbolFactors{1})).AdjClose(2:datalength)./yDataIndex.(genvarname(SymbolFactors{1})).AdjClose(1:datalength-1));

    stdret = (retIndex-mean(retIndex))/std(retIndex);

    initial = [0.0, 1];

    [empCumDis,xi_empCumDis] = ksdensity(stdret(:),'function','cdf');
    [empFreDis,xi] = ksdensity(stdret(:));
    [xlsq, resnorm] = lsqcurvefit( @(x,y)cdf_stdPSM(y, x), initial, xi_empCumDis, empCumDis, [ -0.5, 0.250001], [0.5, Inf]);

    
    pdfPSM = pdf_stdPSM( xi, xlsq );  
    cdfPSM = cdf_stdPSM( xi, xlsq );  

    figure(1),plot(xi, empFreDis, 'k', xi, pdfPSM, 'r', xi, pdfNormal, 'b') 
    figure(2),plot(xi, log(empFreDis), 'k', xi, log(pdfPSM), 'r', xi, log(pdfNormal), 'b') 
    figure(3),plot(xi_empCumDis, empCumDis, 'k', xi_empCumDis, cdfPSM, 'b') 

    [h,pPSM,ksstatPSM,cv] = kstest(stdret(:),[xi_empCumDis(:),cdfPSM(:)]);
    [adPSM, ad2PSM] = ad2test_FFT( stdret(:),[xi_empCumDis(:),cdfPSM(:)] );
    ad2PSM=ad2PSM*length(stdret);
    pAD2PSM=1-AD(length(stdret),ad2PSM);
    
    cdfG = cdf('norm', xi_empCumDis, 0, 1);
    [h,pG,ksstatG,cv] = kstest(stdret(:),[xi_empCumDis(:),cdfG(:)]);
    [adG, ad2G] = ad2test_FFT( stdret(:),[xi_empCumDis(:),cdfG(:)] );
    ad2G=ad2G*length(stdret);
    pAD2G=1-AD(length(stdret),ad2G);
    
    figure(4)
    plot(-4:0.1:4, pdf_stdPSM( -4:0.1:4, xlsq )./pdf('norm', -4:0.1:4, 0, 1));
    figure(5)
    plot(-4:0.1:4, cdf_stdPSM( -4:0.1:4, xlsq )./cdf('norm', -4:0.1:4, 0, 1));
return
