%% Math 315 Lab 6
% The following lab explores the composite trapezoidal rule that is used to
% estimate the integral of different functions.

%% The Code
% The following code computes the composite trapezoidal rule
% approximation for n intervals until the error is less than a tolerance
% level of 1000*epsilon. We estimate the error with the relative difference
% between the current and next approximation of N and 2N intervals
% respectively by assuming $f"(\eta _{2N}) = f"(\eta _{N})$. To compute the
% composite trapezoidal rule approximations for N intervals, we compute the function values
% at equally spaced points. Then we recursively compute the approximation
% for 2N intervals by reusing the values computed in the previous
% approximation and then computing the missing values that are in between
% the previous intervals. To solve the slow convergence issue and when $I_{2N} = I_N$ for small N, we implement
% a maximum and minimum number of times, respectively, that the program
% will loop through and compute the next approximation using 2N intervals.
% During each loop, we also compute and store the error estimate. Using
% these values, we can plot the error as the number of intervals N
% increases.

disp(fileread("composite_trapezoidal_rel.m"));

%% Integral A
% $\int^{2\pi}_{0} e^{sinx} ~dx$
%%
% The composite trapezoidal rule approximations of this integral starts
% near epsilon then jumps up. The start of the graph is likely because
% $I_{2N} = I_N$ for N=2. Afterwards, the error rapidly decreases to
% converge near the tolerance of 1000*epsilon. The rate of convergence is
% $O(N^{-17})$.

close all;
f = @(x) exp(sin(x));
a = 0;
b = 2 * pi;
[i, errs] = composite_trapezoidal_rel(f, a, b);
p = polyfit(log(errs(2:end, 1)), log(errs(2:end, 2)), 1);
x = linspace(min(errs(2:end, 1)), max(errs(2:end, 1)));
y = exp(p(2)) * x .^ p(1);
fprintf('y = %g * x^%g\n', exp(p(2)), p(1));
figure1 = figure();
loglog(errs(:, 1), errs(:, 2), '-o', x, y);
xlabel('N');
ylabel('Relative Error');
legend('Rel Errs', 'Convergence Rate Line');
title('Convergence Graph for Integral A');


%% Integral B
% $\int^{\pi}_{0} e^{sinx} ~dx$
%%
% The composite trapezoidal rule approximations of this integral converge
% in a nearly linear line from the N=4 until the tolerance level. This line
% has a slope of ~-2. Since this line is plotted in log scale, the
% convergence rate of this integral is $O(N^{-2})$

f = @(x) exp(sin(x));
a = 0;
b = pi;
[i, errs] = composite_trapezoidal_rel(f, a, b);
p = polyfit(log(errs(2:end, 1)), log(errs(2:end, 2)), 1);
x = linspace(min(errs(2:end, 1)), max(errs(2:end, 1)));
y = exp(p(2)) * x .^ p(1);
fprintf('y = %g * x^%g\n', exp(p(2)), p(1));
figure2 = figure();
loglog(errs(:, 1), errs(:, 2), '-o', x, y);
xlabel('N');
ylabel('Relative Error');
legend('Rel Errs', 'Convergence Rate Line');
title('Convergence Graph for Integral B');

%% Integral C
% $\int^{2\pi}_{0} cosx ~dx$
%%
% Since this integral = 0, I will use the absolute error instead of
% relative errors. The composite trapezoidal rule approximations of this
% integral converges very rapidly to the tolerance level from N=2 to N=4.
% The rate of convergence is $O(N^{-55})$. The points after are zeros and
% round off errors around epsilon. 

f = @(x) cos(x);
a = 0;
b = 2 * pi;
[i, errs] = composite_trapezoidal_abs(f, a, b);
p = polyfit(log(errs(1:2, 1)), log(errs(1:2, 2)), 1);
x = linspace(min(errs(1:2, 1)), max(errs(1:2, 1)));
y = exp(p(2)) * x .^ p(1);
fprintf('y = %g * x^%g\n', exp(p(2)), p(1));
figure3 = figure();
loglog(errs(:, 1), errs(:, 2), '-o', x, y);
xlabel('N');
ylabel('Relative Error');
legend('Rel Errs', 'Convergence Rate Line');
title('Convergence Graph for Integral C');

%% Integral D
% $\int^{2\pi}_{0} (1+cosx) ~dx$
%%
% The composite trapezoidal rule approximations of this integral converges
% very rapidly to the tolerance level from N=2 to N=4. The rate of
% convergence is $O(N^{-58})$. The points after are zeros and round off
% errors around epsilon. 

f = @(x) 1 + cos(x);
a = 0;
b = 2 * pi;
[i, errs] = composite_trapezoidal_rel(f, a, b);
p = polyfit(log(errs(1:2, 1)), log(errs(1:2, 2)), 1);
x = linspace(min(errs(1:2, 1)), max(errs(1:2, 1)));
y = exp(p(2)) * x .^ p(1);
fprintf('y = %g * x^%g\n', exp(p(2)), p(1));
figure4 = figure();
loglog(errs(:, 1), errs(:, 2), '-o', x, y);
xlabel('N');
ylabel('Relative Error');
legend('Rel Errs', 'Convergence Rate Line');
title('Convergence Graph for Integral D');

%% Integral E
% $\int^{2\pi}_{0} |cosx| ~dx$
%%
% The composite trapezoidal rule approximations of this integral starts
% near epsilon then jumps up. The start of the graph is likely because
% $I_{2N} = I_N$ for N=2. Afterwards, the errors follow a nearly linear
% line until convergence at the tolerance level. This line has a slope of
% ~-2. Since this line is plotted in log scale, the convergence rate of
% this integral is $O(N^{-2})$ 

f = @(x) abs(cos(x));
a = 0;
b = 2 * pi;
[i, errs] = composite_trapezoidal_rel(f, a, b);
p = polyfit(log(errs(3:end, 1)), log(errs(3:end, 2)), 1);
x = linspace(min(errs(3:end, 1)), max(errs(3:end, 1)));
y = exp(p(2)) * x .^ p(1);
fprintf('y = %g * x^%g\n', exp(p(2)), p(1));
figure5 = figure();
loglog(errs(:, 1), errs(:, 2), '-o', x, y);
xlabel('N');
ylabel('Relative Error');
legend('Rel Errs', 'Convergence Rate Line');
title('Convergence Graph for Integral E');

%% Integral F
% $\int^1_0 x ~dx$
%%
% The composite trapezoidal rule approximations of this integral are exact
% since the integral geometrically is a trapezoid and thus converge
% instantly. The errors are not round-off errors but rather small numbers
% added to display a graph in log scale. This indicates that the degree of
% precision is at least 1. 

f = @(x) x;
a = 0;
b = 1;
[i, errs] = composite_trapezoidal_rel(f, a, b);
p = polyfit(log(errs(2:end, 1)), log(errs(2:end, 2)), 1);
x = linspace(min(errs(2:end, 1)), max(errs(2:end, 1)));
y = exp(p(2)) * x .^ p(1);
fprintf('y = %g * x^%g\n', exp(p(2)), p(1));
figure6 = figure();
loglog(errs(:, 1), errs(:, 2), '-o', x, y);
xlabel('N');
ylabel('Relative Error');
legend('Rel Errs', 'Convergence Rate Line');
title('Convergence Graph for Integral F');

%% Integral G
% $\int^1_0 x^2 ~dx$
%%
% The composite trapezoidal rule approximation of this integral converges
% to the tolerance level in a nearly linear line in log scale. This line
% has a slope of ~-2. Since this line is plotted in log scale, the
% convergence rate of this integral is $O(N^{-2})$. This indicates that the
% degree of preicion is less than 2.

f = @(x) x.^2;
a = 0;
b = 1;
[i, errs] = composite_trapezoidal_rel(f, a, b);
p = polyfit(log(errs(2:end, 1)), log(errs(2:end, 2)), 1);
x = linspace(min(errs(2:end, 1)), max(errs(2:end, 1)));
y = exp(p(2)) * x .^ p(1);
fprintf('y = %g * x^%g\n', exp(p(2)), p(1));
figure7 = figure();
loglog(errs(:, 1), errs(:, 2), '-o', x, y);
xlabel('N');
ylabel('Relative Error');
legend('Rel Errs', 'Convergence Rate Line');
title('Convergence Graph for Integral G');

%% Integral H
% $\int^1_0 x^3 ~dx$
%%
% The composite trapezoidal rule approximation of this integral converges
% to the tolerance level in a nearly linear line in log scale. This line
% has a slope of ~-2. Since this line is plotted in log scale, the
% convergence rate of this integral is $O(N^{-2})$ 

f = @(x) x.^3;
a = 0;
b = 1;
[i, errs] = composite_trapezoidal_rel(f, a, b);
p = polyfit(log(errs(2:end, 1)), log(errs(2:end, 2)), 1);
x = linspace(min(errs(2:end, 1)), max(errs(2:end, 1)));
y = exp(p(2)) * x .^ p(1);
fprintf('y = %g * x^%g\n', exp(p(2)), p(1));
figure8 = figure();
loglog(errs(:, 1), errs(:, 2), '-o', x, y);
xlabel('N');
ylabel('Relative Error');
legend('Rel Errs', 'Convergence Rate Line');
title('Convergence Graph for Integral H');

%% Convergence Rates
% From the convergence rates of integrals, F, G, and H, we can determine
% that the degree of precision for the composite trapezoidal rule is 1,
% since the integral for f(x) = x converged instantly but not the integrals
% for f(x) = x^2 and f(x) = x^3. From this, I conjecture that the rate of
% convergence in terms of the number of sub-intervals is $O(N^{-(DOP +
% 1)})$. As the degree of precision is 1, the 
% convergence rate of the integral should be $O(N^{-2})$, which is shown in
% integrals B, E, G, and H. To test this, I will approximate the integral
% of $f(x) = 5x^7 - 8x^3$ on the interval 0 to 1. In the graph shown below,
% the convergence rate is $O(N^{-2})$. 

f = @(x) 5*x.^7-8*x.^3;
a = 0;
b = 1;
[i, errs] = composite_trapezoidal_rel(f, a, b);
p = polyfit(log(errs(2:end, 1)), log(errs(2:end, 2)), 1);
x = linspace(min(errs(2:end, 1)), max(errs(2:end, 1)));
y = exp(p(2)) * x .^ p(1);
fprintf('y = %g * x^%g\n', exp(p(2)), p(1));
figure9 = figure();
loglog(errs(:, 1), errs(:, 2), '-o', x, y);
xlabel('N');
ylabel('Relative Error');
legend('Rel Errs', 'Convergence Rate Line');
title('Convergence Graph for Test Integral 1');

%%
% I also conjecture that integrals over exactly one period of a smooth and 
% periodic function have superconvergence, meaning that it will converge
% much quicker than expected. This is shown in integrals A, C, and D. To
% test this, I will approximate the integral of $f(x) = cosx + sinx + 1$ on
% the interval $\frac{\pi}{4}$ to $\frac{9\pi}{4}$. The convergence rate is
% $O(N^{-52})$ 

f = @(x) sin(x) + cos(x) + 1;
a = pi / 4;
b = 9 * pi / 4;
[i, errs] = composite_trapezoidal_rel(f, a, b);
p = polyfit(log(errs(1:2, 1)), log(errs(1:2, 2)), 1);
x = linspace(min(errs(1:2, 1)), max(errs(1:2, 1)));
y = exp(p(2)) * x .^ p(1);
fprintf('y = %g * x^%g\n', exp(p(2)), p(1));
figure10 = figure();
loglog(errs(:, 1), errs(:, 2), '-o', x, y);
xlabel('N');
ylabel('Relative Error');
legend('Rel Errs', 'Convergence Rate Line');
title('Convergence Graph for Test Integral 2');