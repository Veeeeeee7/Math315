%% n2deriv
% This function approximates the second derivative of a function f using the
% central difference quotient. This approximation is derived by taking the
% taylor series approximation of f and truncating the series after the
% term. The truncation error is of the order $O(h^2)$ and will decrease as
% h decreases.

function dydx2 = n2deriv(f, c, h)
  dydx2 = (f(c-h) - 2.*f(c) + f(c+h)) ./ (h.^2);
end