%% n1deriv 
% This function approximates the first derivative of a function f using the
% central difference quotient. This approximation is derived by taking the
% taylor series approximation of f and truncating the series after the
% term. The truncation error is of the order $O(h^2)$ and will decrease as
% h decreases.

function dydx = n1deriv(f, c, h)
  dydx = (f(c+h) - f(c-h)) ./ (2.*h);
end