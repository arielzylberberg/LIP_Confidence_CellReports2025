function yj = jitter(y, e)

if nargin==1 || isempty(e)
    e = 10^-8;
end

yj = y+rand(size(y))*e;