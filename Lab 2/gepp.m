function x = gepp(A, b)
    %
    % Solves Ax=b using Gaussian elimination with partial pivoting
    % by rows for size. Initializations:
    %
    n = length(b); x = zeros(n,1);
    tol = eps*max(A(:));
    %
    % Loop for stages k = 1, 2, ..., n-1
    %
    for k = 1:n-1
    % Search for pivot entry:
        [piv, psub] = max(abs(A(k:n,k)));
        p = psub + k - 1;

        % Exchange current row, k, with pivot row, p:
        A([k,p],k:n) = A([p,k],k:n);
        b([k,p]) = b([p,k]);

        % Check to see if A is singular:
        if abs(A(k,k)) < tol
            error('Linear system appears to be singular')
        end
    
        % Perform the elimination step - row-oriented:
        A(k+1:n,k) = A(k+1:n,k) / A(k,k);
        for i = k+1:n
            A(i,k+1:n) = A(i,k+1:n) - A(i,k)*A(k,k+1:n);
            b(i) = b(i) - A(i,k)*b(k);
        end
    end
    
    % Check to see if A is singular:
    if abs(A(n,n)) < tol
        error('Linear system appears to be singular')
    end

    % Solve the upper triangular system by row-oriented backward substitution:
    for i = n:-1:1
        x(i) = (b(i) - A(i,i+1:n)*x(i+1:n)) / A(i,i);
    end