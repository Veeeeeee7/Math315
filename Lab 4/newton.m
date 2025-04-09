function [y, b] = newton(x, f, xx)
    n = length(x);
    b = zeros(1, n);
    b(1) = f(1);
    for k = 2:n
        num = f(k) - b(1);
        for j = 2:k-1
            multiplier = (x(k) - x(1));
            for i = 2:j-1
                multiplier = multiplier .* (x(k) - x(i));
            end
            
            num = num - b(j) .* multiplier;
        end

        den = 1;
        for j = 1:k-1
            den = den .* (x(k) - x(j));
        end

        b(k) = num ./ den;
    end

    y = zeros(1, length(xx)) + b(n);
    for i = n-1:-1:1
        % disp(b(i))
        y = b(i) + (xx - x(i)) .* y;
    end
end