function x = col_backward(A, b)
    % Solve Ax=b, where A is an n-by-n upper triangular matrix,
    % using col-oriented backward substitution.
    if any(diag(A) == 0)
        error('Input matrix is singular')
    end
    n = length(b);
    x = zeros(n,1);
    
    for j = n:-1:1
        x(j) = b(j) / A(j, j);
        for i = 1: j-1
            b(i) = b(i) - A(i, j) * x(j);
        end
    end