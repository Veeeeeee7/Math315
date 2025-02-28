%% Logical indexing  piecewise function example
% L01_pwindex(x) computes $1 + \sin 3x$ if $x \le 3$ and $x-5$ if $x > 3$
function y = L01_pwindex(x)
s    = x <= 3;           % creates a logical index array
y    = zeros(size(x));   % preallocate y
y(s) = 1 + sin(3.*x(s)); % store y for x <= 3
s    = not(s);           % equivalent to x > 3
y(s) = x(s) - 5;         % store y for x > 3
end
%%%
% To compare with if-then/arrayfun:
%
%  xx = linspace(0,5,10000);
%  tic; yy=L01_pwindex(xx); toc
%  tic; yy=arrayfun(@(x)L01_pwif(x), xx); toc
%
% For piecewise definitions with many conditions, this approach ends up
% testing all |x| for each condition.  There is no short-circuiting as can
% happen with nested if-else blocks.
