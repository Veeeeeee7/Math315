function [rel_forward_err_inv, rel_backward_err_inv, rel_forward_err_gepp, rel_backward_err_gepp] = rel_err_u(A, U, S, V, Z, k)
    x2 = Z * U(:,k);
    b2 = U * (S * (V' * x2));
    x2z = Z * b2;
    rel_forward_err_inv = norm(x2z - x2) / norm(x2);
    rel_backward_err_inv = norm(A * x2z - b2) / (norm(A) * norm(x2z) + norm(b2));
    x2g = linsolve(A, b2);
    rel_forward_err_gepp = norm(x2g - x2) / norm(x2);
    rel_backward_err_gepp = norm(A * x2g - b2) / (norm(A) * norm(x2g) + norm(b2));

