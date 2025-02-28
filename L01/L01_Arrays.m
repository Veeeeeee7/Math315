%% Example coding
%
%
%% Initialization
%
% Use Run Section on each and examine the results
%
f  = @(x)sin(x);             %
-cos(2)+cos(-1)              % the integral of f(x) from x=-1 to x=2
%% Example 1A
%    Element oriented approach
x  = linspace(-1,2,1001);    % subdivide interval [-1,2]
dx = zeros(1,length(x)-1);   % initialize difference array
for i = (1:length(x)-1)      %
    dx(i) = x(i+1) - x(i);   % store delta x
end                          %
y = zeros(1,length(dx));     % initialize array y
for i = (1:length(dx))       %
    y(i) = f(x(i));          % map f onto array x
end                          %
a1 = 0;                      % initialize accumulator
for i = (1:length(dx))       %
    a1 = a1 + y(i) .* dx(i); % sum y times dx
end                          %
a1                           % output a1
%% Example 1B
%    Array oriented approach
x  = linspace(-1,2,1001);    %
dx = x(2:end)-x(1:end-1);    % equiv. to lines 10-13
y  = f(x(1:end-1));          % equiv. to lines 14-17
a2 = y*dx'                   % equiv. to lines 18-22
%% Example 1C
%    Minor variation on Example 1B
x  = linspace(-1,2,1001);    %
dx = diff(x);                % built-in func. for differences
y  = f(x(1:end-1));          %
a3 = y*dx'                   %
