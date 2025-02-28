%% Math 315: Lab 3
% The following lab uses MATLAB functions to compute norms and condition
% numbers to analyze the difference in error between GEPP and using the
% inverse of a matrix. These errors are then used to check the validity of
% the bounds derived in the problem set.

%% Procedure 1
% In the first procedure, I will randomly generate 10 random matrices A
% with a random vector b and solve Ax=b. Using MatLab's svd function, I
% will create two random orthogonal matrices U and V. I will then discard the
% singular values and generate them using a log scale. Then, I will compute
% the matrix A = U * S * V'. I will then be able to calculate the accurate
% inverse of A, with only double precision round-off errors. Using this 
% accurate inverse, I will be able to calculate the accurate solution x, 
% which can then be used to compare the GEPP and inverse matrix solutions.
% 

n = 256;
sigma_1 = 10^4;
sigma_n = 10^-4;

rel_backward_err_inv = zeros(10, 1);
rel_forward_err_inv  = zeros(10, 1);
rel_backward_err_gepp = zeros(10, 1);
rel_forward_err_gepp  = zeros(10, 1);

saved_A = zeros(256,256,10);

for i = 1:10
    % Generate random matrix A, where the accurate inverse of A can be
    % accurately computed
    [U, ~, V] = svd(randn(n));
    svalues = logspace(log10(sigma_1), log10(sigma_n), n);
    S = diag(svalues);
    invS = diag(svalues.^-1);
    A = U * S * V';
    AccurateInv = V * invS * U';
    Z = inv(A);
    
    % Save A for future use
    saved_A(:,:,i) = A;

    % Generate random vector b and calculate forward/backward error for
    % solution x via inverse of A and GEPP
    b = randn(n, 1);
    x = V * (invS * (U' * b));
    xz = Z * b;
    xg = linsolve(A, b);
    rel_backward_err_inv(i) = norm(A * xz - b) / (norm(A) * norm(xz) + norm(b));
    rel_forward_err_inv(i)  = norm(xz - x) / norm(x);
    rel_backward_err_gepp(i) = norm(A * xg - b) / (norm(A) * norm(xg) + norm(b));
    rel_forward_err_gepp(i)  = norm(xg - x) / norm(x);
end

random_A_errs_inv = table(rel_backward_err_inv, rel_forward_err_inv, ...
    'VariableNames', {'Backward_Err_via_Z', 'Forward_Err_via_Z'});
disp(['Table of Relative Backward and Forward Error of Solutions of 10 Random ' ...
    'Matrices A and vectors b via Inverse Z']);
disp(random_A_errs_inv);

random_A_errs_gepp = table( rel_backward_err_gepp, rel_forward_err_gepp, ...
    'VariableNames', {'Backward_Err_via_GEPP', 'Forward_Err_via_GEPP'});
disp(['Table of Relative Backward and Forward Error of Solutions of 10 Random ' ...
    'Matrices A and vectors b via GEPP']);
disp(random_A_errs_gepp);

random_A_errs_comparison = table(abs(rel_backward_err_inv) ./ abs(rel_backward_err_gepp), ...
    abs(rel_forward_err_inv) ./ abs(rel_forward_err_gepp), ...
    'VariableNames', {'Ratio of Backward Err from Inverse to GEPP', ['Ratio ' ...
    'Forward Err from Inverse to GEPP']});
disp(['Table of Ratios of Relative Backward and Forward Error from Inverse Z ' ...
    'to GEPP for Solutions of 10 Random Matrices A and vectors b. (Ratio = |err ' ...
    'via inverse Z| - |err via GEPP|)']);
disp(random_A_errs_comparison);

%%
% In the tables above, it's evident that solving systems using GEPP generally 
% has lower relative backward error than using the inverse Z as all the values 
% in the table of ratios between relative backward error of inverse Z to 
% that of GEPP are greater than 1. However, there isn't a massive difference
% as the ratio isn't in the orders of magnitude of 10. Therefore, using the 
% inverse Z to solve a linear system will be generally as accurate as GEPP, 
% having around 1 decimal place of accuracy less than GEPP. Looking at the 
% errors themself, it seems like they are all around the magnitude of epsilon 
% meaning that the relative backward error seems to be bounded by epsilon 
% time some small constant. This suggests that the O(machine epsilon) bound 
% for backwards error of GEPP is accurate. 

%%
% To test if the bound is accurate for backwards error of inverse Z, I will 
% calculate the bound using the last matrix calculated and
% compare with the largest relative backward error for inverse Z. As shown
% below the bound for relative backward error for inverse Z is also quite
% accurate.

norm_A_inv = norm(AccurateInv);
norm_xz = norm(xz);
norm_b = norm(b);

fprintf('Bound = O((||A^-1|| * ||b|| * 2*10^-16) / (||xz||) ~= %e\n', ...
    norm_A_inv .* norm_b .* eps() ./ norm_xz);
fprintf('Largest actual error ~= %e\n', max(abs(rel_backward_err_inv)));

%%
% For relative forward error, the ratio is almost nearly one for every entry. 
% This follows from equations (1.5) and (1.6) that show the bounds of relative 
% forward error for solutions via GEPP and inverse Z are equal if Z is a 
% sufficiently accurate inverse.

%%
% I will use the calculated condition numbers below for the different matrices 
% A to check if the bounds described in (1.5) and (1.6) are reasonable. Since 
% they are all approximately 10^8 as shown below and the error bounds are
% approximately 10^8 * 2 * 10^-16 = 2 * 10^-8, which is just greater than the
% calculated relative errors. Hence, the bounds described in equations (1.5) 
% and (1.6) are quite accurate.

for i=1:10
    c = cond(saved_A(:,:,i));
    fprintf('Matrix %d condition number: %e\n', i, c);
end

%% Procedure 2
% In the second procedure, I will randomly generate a single matrix A using
% MatLab's svd function in the same way as procedure 1 to ensure an accurate 
% inverse can be calculated. I will then solve the equation Ax=uj where uj
% equals each of the left singular vectors of A. Using the solution from
% the accurate inverse, I will then calculate the relative backward and 
% forward error of these solutions.

% Generate random matrix A, where the accurate inverse of A can be
% accurately computed
n = 256;
sigma_1 = 10^4;
sigma_n = 10^-4;
[U2, ~, V2] = svd(randn(n));
svalues2 = logspace(log10(sigma_1), log10(sigma_n), n);
S2 = diag(svalues2);
invS2 = diag(svalues2.^-1);
A2 = U2 * S2 * V2';
AccurateInv2 = V2 * invS2 * U2';
Z2 = inv(A2);

%%
% For each left singular vector in U, I will solve Ax=uj and calculate 
% relative backward and forward error using the rel_err_u() function 
% displayed below.

disp(fileread("rel_err_u.m"));

%%
rel_backward_err_inv2 = zeros(256, 1);
rel_forward_err_inv2  = zeros(256, 1);
rel_backward_err_gepp2 = zeros(256, 1);
rel_forward_err_gepp2  = zeros(256, 1);
for i = 1:256
    [rel_forward_err_inv2(i), rel_backward_err_inv2(i), rel_forward_err_gepp2(i), ...
        rel_backward_err_gepp2(i)] = rel_err_u(A2, U2, S2, V2, Z2, i);
end

figure();
loglog(svalues2, rel_backward_err_inv2, svalues2, rel_backward_err_gepp2);
xlabel('Singular Value \sigma_j Corresponding to u_j');
ylabel('Relative Backward Error of Ax=u_j');
legend('Inverse Z', 'GEPP');
title('Relative Backward Error vs Singular Value \sigma_j')

figure();
loglog(svalues2, rel_forward_err_inv2, svalues2, rel_forward_err_gepp2);
xlabel('Singular Value \sigma_j Corresponding to u_j');
ylabel('Relative Forward Error of Ax=u_j');
legend('Inverse Z', 'GEPP');
title('Relative Forward Error vs Singular Value \sigma_j')

%%
% From the first graph, we can see a major difference in relative backward 
% error between using inverse Z and GEPP. Nearly 6 digits of accuracy were 
% lost when the vector corresponding to the largest singular value was used. 
% This is because the bound for relative backward error for inverse Z can be 
% simplified to O(corresponding singular value / minimum singular value) as 
% shown in problem #5 in the attached PDF. This formula suggests that as the 
% singular values corresponding to the left singular vectors used for b
% increase, the bound on relative backward error will also increase. This
% is reflected in the graph as the relative errors follow a pretty
% linearly increasing line as the corresponding singular value increases.

%%
% To check the accuracy of this bound, I will calculate the ratio of the
% largest and smallest singular value, then multiplied by machine epsilon.
% Comparing this to the actual largest backward error from inverse Z will 
% show that this bound is accurate.

fprintf('Bound = %e\n', (sigma_1 ./ sigma_n) .* eps());
fprintf('Largest actual error = %e\n', max(abs(rel_backward_err_inv2)));

%%
% Also from the first graph, the relative backward error of GEPP is around 
% machine epsilon as predicted by the bound of O(machine epsilon).
% Therefore this error bound for GEPP is accurate.

%% 
% Lastly, the second graph showcasing relative forward error doesn't have a
% significant difference between the error from using inverse Z and using 
% GEPP. The errors are generally of one magnitude of 10 or less, meaning 
% that at most one digit of accuracy was lost when solving using inverse Z 
% as compared to GEPP. From the condition number of A calculated below, we
% can approximate the bound to be around 10^8 * 2 * 10^-16 = 2 * 10^-8. 
% Hence, the error bound for forward error is also accurate.

fprintf('Condition number of A: %e\n', cond(A2));
fprintf('Largest actual error of inverse Z = %e\n', max(abs(rel_forward_err_inv2)));
fprintf('Largest actual error of GEPP = %e\n', max(abs(rel_forward_err_gepp2)));