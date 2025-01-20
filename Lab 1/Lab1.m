x = linspace(0,2*pi,100);
h=0.001;
approx = nderiv(x,h);
actual = sincP(x);
err = actual - approx;
disp(err)