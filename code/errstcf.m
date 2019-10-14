function f = errstcf(x, samplelogcf, wk, functionhandle)
    
    theory = log(functionhandle(wk, x));
    f = [real(theory) - real(samplelogcf); imag(theory) - imag(samplelogcf)];
    f = f*100;
return