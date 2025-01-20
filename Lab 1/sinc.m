function f = sinc(x)
f = sin(x) ./ x;
f(x == 0) = 1;
end