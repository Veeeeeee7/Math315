function f = nderiv(x,h)
f = (sinc(x+h) - sinc(x)) ./ h
end