function f = sincP(x)
solution = zeros(1,length(x));
for i = 1:length(x)
    if x(i)~=0
        solution(i) = (x(i) * cos(x(i)) - sin(x(i))) / (x(i)^2);
    end
f = solution;
end
     