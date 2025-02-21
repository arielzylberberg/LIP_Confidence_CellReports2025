function [c2,gof,vals,fallo,x] = fit_PRR_function(t,y,varargin)
% [c2,gof,vals] = fit_PRR_function(t,y)
% t: time
% pars: mu,sigma,alfa,c,d

% parameters: alfa,c,d,mu,sigma
x0 = [0.2  1 0.5 .2  .2];
for i=1:length(varargin)
    if isequal(varargin{i},'x0')
        x0 = varargin{i+1};
    end
end

if size(t,2) > size(t,1)
    t = t';
end

options = optimset('MaxIter',3000,'MaxFunEvals',1500,'TolFun',1e-8);
%options = optimset('MaxIter',3000,'MaxFunEvals',15000,'TolFun',1e-8);
lb = [-1 0 0 0 0];
ub = ones(1,5) * 100;
lb(2) = 1;
ub(2) = 1;

[x,resnorm,residual,exitflag] = lsqnonlin(@myfun,x0,lb,ub,options);


vals.t = t;
ymodel = prr_function(x(1),x(2),x(3),x(4),x(5));
vals.fit = ymodel;
c2 = [];
gof = RSquared(y,ymodel);
fallo = not(exitflag==1);

    function F = myfun(x)
        ymodel = prr_function(x(1),x(2),x(3),x(4),x(5));
        F = ymodel - y;
    end

    function ymodel = prr_function(alfa,c,d,mu,sigma)
        
        % old:
        % ymodel = d*exp(mu*alfa + 0.5*sigma^2*alfa^2 - ...
        %     alfa*t).*normcdf(t,mu+sigma^2*alfa,sigma) + c*normcdf(t,mu,sigma);


        ymodel = d * (c * exp(mu*alfa + 0.5*sigma^2*alfa^2 - ...
            alfa*t).*normcdf(t,mu+sigma^2*alfa,sigma)) + (1-c)*normcdf(t,mu,sigma);

    end

end