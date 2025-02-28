%% Code examples

%% plot()
% plot, line specification

xx = linspace(0, 2, 100);
yy = exp(xx); % array evaluation of e^x
fig1 = figure;
plot(xx, yy, '--') % dashed

%% semilogy()
% "log-plot", modify properties

xx = linspace(0, 2, 20);
yy = exp(2 .* xx); % array evaluation of e^(2x)
fig2 = figure;
myplot = semilogy(xx, yy)
hold on % after log plot
myplot.LineStyle = ":"; % dotted
myplot.LineWidth = 2; % dotted
myplot.Color = "magenta";
myplot.Marker = ".";
myplot.MarkerSize = 20;
hold off

%% loglog()
% "log-log-plot", markers, legend

xx = logspace(1, 5, 15);
y2 = xx.^2; % array evaluation of x^2
y3 = xx.^3; % array evaluation of x^3
fig3 = figure;
loglog(xx, y2, '*', xx, y3, 'o') % point markers '*', 'o'
hold on % after log plot
legend('x^2','x^3','Location','northwest')
hold off

%% Logical indexing

xx = linspace(0,2.*pi,16);
ss = sin(xx);
cc = cos(xx);
ctX = cc/ss;   % A scalar? What does this do??? (WRONG)
ctg = cc./ss;  % cotangent(x) (RIGHT)
Defined = isfinite(ctg); % logical array
ctg(Defined); % values in array ctg for which Defined is 1 (true)
max(abs(ctg(Defined))) % max defined/finite value
max(abs(ctg))          % max of all values (incl. inf)