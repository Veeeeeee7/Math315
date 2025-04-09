%=======================================
% DCT-1 : Discrete cosine transform
%---------------------------------------
% Code adapted from
% P. Moin, Fundamentals of Engineering Numerical Analysis (2010)
%
% Computes the same result as exemplified in Chebseries.m
% but uses the built-in Fast Fourier Transform fft()
%
% The following are equivalent:
%    c = DCT1(exp(cos((0:16)./16.*pi))')
%    c = Chebseries(@(x) exp(x), 16)
%
%=======================================
function y=dct1(f_vals)

N      =length(f_vals);
y      =[f_vals;flipud(f_vals(2:N-1,:))];
y_ft   =fft(y);
y      =real(y_ft(1:N,:)/(N-1));
y(1,:) = y(1,:)/2;
y(N,:) = y(N,:)/2;
