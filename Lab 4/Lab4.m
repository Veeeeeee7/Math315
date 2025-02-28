warning('off', 'MATLAB:nearlySingularMatrix');
n=100;
x=linspace(-1,1,n);
y=exp(x);
v=vander(x);
p=v\y';
t=linspace(-1,1,1025.*n);
err=max(abs(polyval(p,t)-exp(t)));
disp(err);

myt = linspace(-0.5, 0.5, 1000);
plot(myt, polyval(p,myt));