function [u,y] = fold_for_logplot(coh,v,type)

if nargin<3
    type = 'data';
end

if nanmax(coh)>1
    coh = coh/100;
end

acoh = abs(coh);
first_positive = 0.032;
zero_val = first_positive/2;


[u,y] = curva_media(v,acoh,[],0);
I = u>first_positive*0.75;
y(~I) = nan;

% add zero and horizontal line

x2 = zero_val;
if isequal(type,'model')
    x3 = 1.2*zero_val;
    x1 = exp(2*log(x2)-log(x3));
    u = [x1; x2; x3; u(:)];
    val = nanmean(v(acoh==0));
    y = [val; val; val; y(:)];
else
    u = [x2; u(:)];
    val = nanmean(v(acoh==0));
    y = [val; y(:)];
end

