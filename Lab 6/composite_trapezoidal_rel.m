function [is, errs] = composite_trapezoidal_rel(f, a, b)
    tol = 1000 * eps;
    min = 2^5;
    max = 10.^8;
    dx = b-a;
    i = [0.5, 0.5] * f([a; b]) .* dx;
    n = 1;
    is = [];
    errs = [];
    while n < min || (abs((i-i2) ./ i2) > tol && n < max)
        dx = dx ./ 2;
        di = sum(f(a+dx:2*dx:b-dx)) * dx;
        i2 = i;
        i = i ./ 2 + di;
        n = 2*n;
        is = [is; n, i];
        errs = [errs; n, abs((i-i2) ./ i2)];
    end
    errs = errs + tol/100000;
end