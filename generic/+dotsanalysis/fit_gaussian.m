function [f1,xx,yy] = fit_gaussian(x,y,startPoints,doplot)

if nargin<3 || isempty(startPoints)
    startPoints = [1 100 0];
end
% gaussEqn = 'a*exp(-((x-b)/c)^2)+d';
% startPoints = [1.5 900 10 0.6];

gaussEqn = 'a*exp(-((x)/c)^2)+d';

if nanmax(x)<1
    norm = 1000;
else
    norm = 1;
end
x = x*norm;

f1 = fit(x,y,gaussEqn,'Start', startPoints);
xx = linspace(min(x),max(x),100);
yy = f1(xx);
xx = xx/norm;
if doplot
    plot(xx,yy);
end