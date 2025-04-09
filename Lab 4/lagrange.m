function [y, l_saved] = lagrange(x, f, xx)
    y = 0;
    n = length(x);
    l_saved = zeros(1, n);
    for k = 1:n
        l_num = 1;
        l_den = 1;
        for j = 1:k-1
            l_num = l_num .* (xx - x(j));
            l_den = l_den .* (x(k) - x(j));
        end
        for j = k+1:n
            l_num = l_num .* (xx - x(j));
            l_den = l_den .* (x(k) - x(j));
        end
    l = l_num ./ l_den;
    if isscalar(l)
        l_saved(k) = l;
    else
        l_saved(k) = l(k);
    end
    y = y + f(k) * l;
    end
end
