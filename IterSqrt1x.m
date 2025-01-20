%% Iteration: x = sqrt(1+x)
%
x = 1;
y = 2;
phi = (1+sqrt(5)) ./2;
while x ~= y
    y = x;
    x = sqrt(1+x);
    err = x-phi;
    err_ratio = err./(y-phi)
end
x