function [ d ] = DividedDifference(x, d)
% DividedDifference computes the Newton divided differences
%   x = list of x input numbers
%   d = list of function values at x (intially)
%   d = current divided differences (computed in place)

n=length(x);
for k=2:n
    for j=n:-1:k
        d(j) = (d(j) - d(j-1)) / (x(j) - x(j-k+1));
    end
end

end

