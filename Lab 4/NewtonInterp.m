function [ p ] = NewtonInterp( x, d, t )
%NewtonInterp interpolates Newton divided differences
%   x = list of interpolation points
%   d = list of divided differences of a function f
%   t = input for f(t)
n=length(x);
p=zeros(1,length(t));
p=p+d(n);
for k=n-1:-1:1
    p = d(k) + (t - x(k)) .* p;
end


end

