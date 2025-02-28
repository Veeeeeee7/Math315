function x = row_backward(A, b)
    % Solve Ax=b, where A is an n-by-n upper triangular matrix,
    % using row-oriented backward substitution.
    n = length(b);
    x = zeros(n, 1);
    if any(diag(A) == 0)
        error('Input matrix is singular')
    end
    for i = n:-1:1
        for j = i+1:n
            b(i) = b(i) - A(i, j) * x(j);
        end
        x(i) = b(i) / A(i, i);
    end
end