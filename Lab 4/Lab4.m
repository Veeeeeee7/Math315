%% Math 315 Lab 4
% The following lab examines the Runge phenomenon for 5 different ways of
% polynomial interpolation: Vandermonde / Power Series, Newton, Newton with
% Divided Differences, Lagrange, and Chebyshev interpolation. I will use
% these different methods of interpolation to approximate the exponential
% function and calculate the maximum absolute error of these approximations
% on the interval [-1, 1].

%% Vandermonde / Power Series Interpolation
% The following code takes in nodes (x, f) and approximation points xx to
% outputs approximated values y by creating an interpolating polynomial
% through Vandermonde / Power Series interpolation.

disp(fileread('vandermonde.m'));

%% Newton Textbook Interpolation
% The following code takes in nodes (x, f) and approximation points xx to
% outputs approximated values y by creating an interpolating polynomial
% through the textbook version of Newton interpolation.

disp(fileread('newton.m'));

%% Newton Divided Differences Interpolation
% The following code takes in nodes (x, f) and approximation points xx to
% outputs approximated values y by creating an interpolating polynomial
% through the divided differences version of Newton interpolation.

disp(fileread('newton_divided_diff.m'));

%% Lagrange Interpolation
% The following code takes in nodes (x, f) and approximation points xx to
% outputs approximated values y by creating an interpolating polynomial
% through Lagrange interpolation.

disp(fileread('lagrange.m'));

%% Chebyshev Interpolation
% The following code takes in nodes (x, f) and approximation points xx to
% outputs approximated values y by creating an interpolating polynomial
% through Chebyshev interpolation.

disp(fileread('chebyshev.m'));

%% Equally Spaced Data Points on the Interval [-1, 1]
% The following code generates 100 equally spaced nodes (x, e^x) from -1 to
% 1. It will then create an interpolating polynomial using the different
% interpolation methods. Then, 102500 sample points are then chose to test
% the interpolating polynomial and the maximum absolute error is
% calculated and put together in the table below.

close all;
warning('off', 'MATLAB:nearlySingularMatrix');
vandermonde_err = zeros(10, 1);
newton_err = zeros(10, 1);
newton_divided_diff_err = zeros(10, 1);
lagrange_err = zeros(10, 1);
chebyshev_err = zeros(10, 1);

for n = 10:10:100
    x = linspace(-1, 1, n)';
    f = exp(x);
    xx = linspace(-1, 1, 1025.*n);
    t = linspace(-1, 1, 1000);
    
    % Vandermonde / Power Series Interpolation
    err = max(abs(vandermonde(x, f, xx) - exp(xx)));
    vandermonde_err(n ./ 10) = err;
    % disp(err);
    % f1 = figure(1);
    % plot(t, vandermonde(x, f, t));
    
    % Newton Textbook Interpolation
    err = max(abs(newton(x, f, xx) - exp(xx)));
    newton_err(n ./ 10) = err;
    % disp(err);
    % f2 = figure(2);
    % plot(t, newton(x, f, t));
    
    % Netwon with Divided Differences Interpolation 
    err = max(abs(newton_divided_diff(x, f, xx) - exp(xx)));
    newton_divided_diff_err(n ./ 10) = err;
    % disp(err);
    % f3 = figure(3);
    % plot(t, newton_divided_diff(x, f, t));
    
    % Lagrange Interpolation
    [y, l] = lagrange(x, f, xx);
    err = max(abs(y - exp(xx)));
    lagrange_err(n ./ 10) = err;
    % disp(err);
    % f4 = figure(4);
    % plot(t, lagrange(x, f, t));
    
    % Chebyshev Interpolation
    [y, T] = chebyshev(x, f, xx);
    err = max(abs(y - exp(xx)));
    chebyshev_err(n ./ 10) = err;
    % disp(err);
    % f5 = figure(5);
    % plot(t, chebyshev(x, f, t));
end

vandermonde_err = categorical(compose('%.7e', round(vandermonde_err, 7, ...
    'significant')));
newton_err = categorical(compose('%.7e', round(newton_err, 7, 'significant')));
newton_divided_diff_err = categorical(compose('%.7e', round(newton_divided_diff_err, ...
    7, 'significant')));
lagrange_err = categorical(compose('%.7e', round(lagrange_err, 7, 'significant')));
chebyshev_err = categorical(compose('%.7e', round(chebyshev_err, 7, 'significant')));

T1 = table(linspace(10,100,10)', vandermonde_err, newton_err, newton_divided_diff_err, ...
    'VariableNames', {'n', 'Vandermonde', 'Newton', 'Newton Divided Difference'});
disp('Table 4.3: Errors for Fitting e^x with Equally Spaced Nodes')
disp(T1);
T2 = table(linspace(10,100,10)', lagrange_err, chebyshev_err, 'VariableNames', ...
    {'n', 'Lagrange', 'Chebyshev'});
disp('Table 4.4: Errors for Fitting e^x with Chebyshev Nodes')
disp(T2);


%% Chebyshev Points on the Interval [-1, 1]
% The following code generates 100 nodes (x, e^x) using the Chebyshev
% method for generating points. This distribution will have more points
% near the ends of the interval and less near the center. Then, using the
% different interpolation methods, an interpolating polynomial is created.
% These interpolating polynomials are then tested with 102500 test points
% and the maxmimum absolute error is calculated. These maximums are then
% put together in the table below.

close all;
warning('off', 'MATLAB:nearlySingularMatrix');
vandermonde_err = zeros(10, 1);
newton_err = zeros(10, 1);
newton_divided_diff_err = zeros(10, 1);
lagrange_err = zeros(10, 1);
chebyshev_err = zeros(10, 1);

for n = 10:10:100
    i = linspace(1, n, n);
    a = -1;
    b = 1;
    x = (((b + a) ./ 2) - ((b - a) ./ 2) .* cos((2 .* i + 1) .* pi ./ ...
        (2 .* n + 2)))';
    
    f = exp(x);
    xx = linspace(-1, 1, 1025.*n);
    t = linspace(-1, 1, 1000);
    
    % Vandermonde / Power Series Interpolation
    err = max(abs(vandermonde(x, f, xx) - exp(xx)));
    vandermonde_err(n ./ 10) = err;
    % disp(err);
    % f1 = figure(1);
    % plot(t, vandermonde(x, f, t));
    
    % Newton Textbook Interpolation
    [y, b] = newton(x, f, xx);
    err = max(abs(y - exp(xx)));
    newton_err(n ./ 10) = err;
    % disp(err);
    % f2 = figure(2);
    % plot(t, newton(x, f, t));
    
    % Netwon with Divided Differences Interpolation 
    [y, d] = newton_divided_diff(x, f, xx);
    err = max(abs(y - exp(xx)));
    newton_divided_diff_err(n ./ 10) = err;
    % disp(err);
    % f3 = figure(3);
    % plot(t, newton_divided_diff(x, f, t));
    
    % Lagrange Interpolation
    [y, l] = lagrange(x, f, xx);
    err = max(abs(y - exp(xx)));
    lagrange_err(n ./ 10) = err;
    % disp(err);
    % f4 = figure(4);
    % plot(t, lagrange(x, f, t));
    
    % Chebyshev Interpolation
    [y, T] = chebyshev(x, f, xx);
    err = max(abs(y - exp(xx)));
    chebyshev_err(n ./ 10) = err;
    % disp(err);
    % f5 = figure(5);
    % plot(t, chebyshev(x, f, t));
end

vandermonde_err = categorical(compose('%.7e', round(vandermonde_err, 7, ...
    'significant')));
newton_err = categorical(compose('%.7e', round(newton_err, 7, 'significant')));
newton_divided_diff_err = categorical(compose('%.7e', round(newton_divided_diff_err, ...
    7, 'significant')));
lagrange_err = categorical(compose('%.7e', round(lagrange_err, 7, 'significant')));
chebyshev_err = categorical(compose('%.7e', round(chebyshev_err, 7, 'significant')));

T1 = table(linspace(10,100,10)', vandermonde_err, newton_err, newton_divided_diff_err, ...
    'VariableNames', {'n', 'Vandermonde', 'Newton', 'Newton Divided Difference'});
disp(T1);
T2 = table(linspace(10,100,10)', lagrange_err, chebyshev_err, 'VariableNames', ...
    {'n', 'Lagrange', 'Chebyshev'});
disp(T2);

%% Textbook Newton vs. Divided Difference Newton Interpolation Errors
% The textbook implementation of Newton interpolation creates the least
% accurate interpolating polynomial of e^x for both equally spaced and
% Chebyshev nodes. The divided difference implementation of Newton
% interpolation creates a much better performing polynomial that
% approximates e^x. It ties for second worst with the Vandermonde / Power
% Series interpolating polynomial for equally spaced nodes, however it is
% worse than the Vandermonde / Power Series interpolating polynomial for
% Chebyshev nodes. This is because the textbook version uses the previously
% calculated polynomial to calculate the coefficient whereas the divided
% difference method recursively uses the previous two divided differences
% to calculate the next one. We can see this by calculating the differences
% between the coefficients from the textbook Newton interpolation and the
% divided differences from the divided differences version of Newton
% interpolation.

x = linspace(-1, 1, 100)';
f = exp(x);
xx = 0.994;

[y, b] = newton(x, f, xx);
[y, d] = newton_divided_diff(x, f, xx);

max_diff = 0;
b_m = 0;
d_m = 0;
for i=1:length(b)
    if abs(b(i) - d(i)) > max_diff
        b_m = b(i);
        d_m = d(i);
    end
end

fprintf(['Difference between polynomial and divided difference Newton ' ...
    'interpolation: %e\n'], abs((b_m - d_m) ./ d_m));

%%
% This error is around 10^11. Assuming the smallest (x-xi) ~= 0.01 = 10^-2,
% the error is around 10^9, which is approximately the difference seen between the
% maximum error of textbook's Newton interpolation and the
% divided-difference Newton interpolation at n=100. 

%%
% The same process can be done with the Chebyshev points to compute the
% maximum difference in error between the textbook implementation and
% divided differences implementation of Newton interpolation.

n = 100;
i = linspace(1, n, n);
a = -1;
b = 1;
x = (((b + a) ./ 2) - ((b - a) ./ 2) .* cos((2 .* i + 1) .* pi ./ ...
    (2 .* n + 2)))';
f = exp(x);
xx = 0.994;

[y, b] = newton(x, f, xx);
[y, d] = newton_divided_diff(x, f, xx);

max_diff = 0;
b_m = 0;
d_m = 0;
for i=1:length(b)
    if abs(b(i) - d(i)) > max_diff
        b_m = b(i);
        d_m = d(i);
    end
end

fprintf(['Difference between polynomial and divided difference Newton ' ...
    'interpolation: %e\n'], abs((b_m - d_m) ./ d_m));

%%
% The result above is unreasonably large compared to the difference between
% the errors of the textbook and divided differences implementation of
% Newton interpolation. This suggests that likely an equally large negative
% difference in coefficient could have occured to cancel ou this large
% positive difference.


%% Vandermonde / Power Series Interpolation Error
% The Vandermonde / Power Series interpolating polynomial ties for the 
% second worst approximation of e^x for equally spaced nodes and is 3rd 
% best for Chebyshev points. It has an error bounded by the condition
% number of the Vandermonde matrix times machine epsilon.

x = linspace(-1, 1, 100)';
i = linspace(1, 100, 100);
a = -1;
b = 1;
x_c = (((b + a) ./ 2) - ((b - a) ./ 2) .* cos((2 .* i + 1) .* pi ./ ...
    (2 .* 100 + 2)))';
fprintf('Error bound for equally spaced points Vandermonde: %e\n', ...
    cond(vander(x)) * eps());
fprintf('Error bound for Chebyshev points Vandermonde: %e\n', ...
    cond(vander(x_c)) * eps());


%%
% For the Chebyshev points, the error is 3 decimal places below the error bound for the
% GEPP. This is because the calculation of the Chebyshev points loses 2
% decimal places of accuracy. On the other hand, for the equally spaced points, the error is above the error
% bound for GEPP. This is likely because of the Runge phenomenon that
% states the error for an interpolating polynomial increases arbitrarily as
% n approaches infinity.

%% Lagrange Interpolation Error
% The Lagrange interpolating polynomial has the second best approximation
% of e^x for both equally spaced nodes and Chebyshev points. We can
% calculate the error bound for Lagrange interpolation by the formula
% below.

x = linspace(-1, 1, 100)';
i = linspace(1, 100, 100);
a = -1;
b = 1;
x_c = (((b + a) ./ 2) - ((b - a) ./ 2) .* cos((2 .* i + 1) .* pi ./ ...
    (2 .* 100 + 2)))';
xx = 0.994;

[y, l] = lagrange(x, f, xx);
[y, l_c] = lagrange(x_c, f, xx);

err_bound = max(abs(l)) * eps();
err_bound_c = max(abs(l_c)) * eps();

fprintf('Error bound for equally spaced points Lagrange Interpolation: %e\n', ...
    err_bound);
fprintf('Error bound for Chebyshev points Lagrange Interpolation: %e\n', ...
    err_bound_c);

%%
% The error bound for equally spaced points is quite accurate whereas the
% error bound for Chebyshev points about 3 decimal places off. This is because 
% the calculation of the Chebyshev points loses 2 decimal places of accuracy.

%% Chebyshev Interpolation Error
% The Chebyshev interpolating polynomial has the best approximation of e^x
% for both equally spaced nodes and Chebyshev points. It's error can be
% bounded by the error of the GEPP subprocess that occurs in the
% coefficient calculations. The error will be bounded by the condition
% number of the matrix T times machine epsilon.

x = linspace(-1, 1, 100)';
i = linspace(1, 100, 100);
a = -1;
b = 1;
x_c = (((b + a) ./ 2) - ((b - a) ./ 2) .* cos((2 .* i + 1) .* pi ./ ...
    (2 .* 100 + 2)))';
xx = 0.994;
[y, T] = chebyshev(x, f, xx);
[T_c,y_c] = chebyshev(x_c, f, xx);
fprintf('Error bound for equally spaced points Chebyshev Interpolation: %e\n', ...
    cond(T) * eps());
fprintf('Error bound for Chebyshev points Chebyshev Interpolation: %e\n', ...
    cond(T_c) * eps());

%%
% The error bound for equally spaced points is quite accurate which suggests 
% that the main source of error was from the GEPP subprocess that is used to 
% calculate the coefficients of the polynomial. However the error bound for
% Chebyshev points is around 3 decimal places off. This is because the 
% calculation of the Chebyshev points loses 2 decimal places of accuracy.