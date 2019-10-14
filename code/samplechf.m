function f = samplechf(u, sample)
    f = zeros(size(u));
    for k = 1:length(u)
        f(k) = mean(exp(1i*u(k)*sample));
    end
return