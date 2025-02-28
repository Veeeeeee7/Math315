function x = row_forward(A, b)
    % Solve Ax=b, where A is an n-by-n lower triangular matrix,
    % using row-oriented forward substitution.
    if any(diag(A) == 0)
        error('Input matrix is singular')
    end
    n = length(b);
    x = zeros(n,1);
    for i = 1:n
        x(i) = (b(i) - A(i,1:i-1)*x(1:i-1)) / A(i,i);
    end