function [ y, d ] = newton_divided_diff( x, f, xx )
    n = length(x);
    d = f;
    for k = 2:n
        for j = n:-1:k
            
            d(j) = (d(j) - d(j-1)) / (x(j) - x(j-k+1));
            % disp(d(j));
        end
    end

    y = zeros(1, length(xx));
    y = y+d(n);
    for k = n-1:-1:1
        y = d(k) + (xx - x(k)) .* y;
    end
end