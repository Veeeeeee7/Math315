function [y] = vandermonde(x, f, xx)
    v = vander(x);
    p = v\f;
    y = polyval(p,xx);
end