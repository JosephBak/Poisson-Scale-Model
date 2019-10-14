SymbolFactors = { '^GSPC' };
yDataIndex = getYahooDailyData(SymbolFactors, '01/01/2013', '12/31/2015', 'mm/dd/yyyy');
datalength = length(yDataIndex.(genvarname(SymbolFactors{1})).AdjClose);
numoffactor = length(SymbolFactors);
retIndex = ...
    log(yDataIndex.(genvarname(SymbolFactors{1})).AdjClose(2:datalength)./yDataIndex.(genvarname(SymbolFactors{1})).AdjClose(1:datalength-1));

    stdret = (retIndex-mean(retIndex))/std(retIndex);

    wk = 0.1:0.1:1;
    wk = wk(:);
    samplelogsymcf = log(samplechf(wk, stdret));

    initial = [0.0, 1];
    opt = optimoptions(@lsqnonlin, 'TolX', 1e-9, 'TolFun', 1e-9, 'MaXFunEvals', 1000);
    [xlsq, resnorm] = lsqnonlin(@(x)errstcf(x, samplelogsymcf, wk, @chf_stdPSM), initial, [-0.5,0.250001], [0.5,inf], opt);
    
    [empCumDis,xi_empCumDis] = ksdensity(stdret(:),'function','cdf');
    [empFreDis,xi] = ksdensity(stdret(:));
    pdfPSM = pdf_stdPSM( xi, xlsq );  
    cdfPSM = cdf_stdPSM( xi, xlsq );  

    pdfNormal = pdf('norm', xi, 0, 1);
    cdfNormal = cdf('norm', xi, 0, 1);

    figure(1),plot(xi, empFreDis, 'k', xi, pdfPSM, 'r', xi, pdfNormal, 'b') 
    figure(2),plot(xi, log(empFreDis), 'k', xi, log(pdfPSM), 'r', xi, log(pdfNormal), 'b') 
    figure(3),plot(xi_empCumDis, empCumDis, 'k', xi_empCumDis, cdfPSM, 'b') 

    [h,pPSM,ksstatPSM,cv] = kstest(stdret(:),[xi_empCumDis(:),cdfPSM(:)]);
    [adPSM, ad2PSM] = ad2test_FFT( stdret(:),[xi_empCumDis(:),cdfPSM(:)] );
    ad2PSM=ad2PSM*length(stdret);
    pAD2PSM=1-AD(length(stdret),ad2PSM);
    
    
    figure(4)
    plot(-4:0.1:4, pdf_stdPSM( -4:0.1:4, xlsq )./pdf('norm', -4:0.1:4, 0, 1));
    figure(5)
    plot(-4:0.1:4, cdf_stdPSM( -4:0.1:4, xlsq )./cdf('norm', -4:0.1:4, 0, 1));
return
