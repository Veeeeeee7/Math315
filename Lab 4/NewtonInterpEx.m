%% Newton interpolation example
%
func = @(x) x .* sin(4.*x);
n = input('Enter n: ');
x=linspace(-1,1,n);
t=linspace(-1,1,8.*n);
y=func(x);
yt = func(t);
d=y;
d=DividedDifference(x,d);
ny=NewtonInterp(x,d,t);
max(abs(ny - yt))
