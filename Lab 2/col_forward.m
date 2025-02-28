function x = col_forward(A, b)
    % Solve Ax=b, where A is an n-by-n lower triangular matrix,
    % using column-oriented forward substitution.
    if any(diag(A) == 0)
        error('Input matrix is singular')
    end
    n = length(b);
    x = zeros(n,1);
    
    for j = 1:n
        x(j) = b(j) / A(j,j);
        b(j+1:n) = b(j+1:n) - A(j+1:n,j)*x(j);
    end