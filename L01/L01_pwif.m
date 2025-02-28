%% If-else piecewise function example
% L01_pwif(x) computes $1 + \sin 3x$ if $x \le 3$ and $x-5$ if $x > 3$
function y = L01_pwif(x)
if x <= 3
    y = 1 + sin(3.*x); % store y for x <= 3
else
    y = x - 5;         % store y for x > 3
end
%%%
% The code |if x <= 3| works only if |x| is a scalar. To apply
% |L01_pwif()| to an array |x|, use |arrayfun()|. To compare with indexing:
%
%  xx = linspace(0,5,10000);
%  tic; yy=L01_pwindex(xx); toc
%  tic; yy=arrayfun(@(x)L01_pwif(x), xx); toc
%