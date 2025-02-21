function y = shuffle_percat(v,cat,nshuffles)
% shuffles the content of vector v, but within 
% category specified by cat. 
% cat: is either a vector, or a concatenation of column vectors
% in which case it unique(rows) is used
% nshuffles: number of vector shuffles to do

if nargin<3
    nshuffles = 1;
end
if isvector(cat)
    c = cat(:);
else
    [~,~,c] = unique(cat,'rows');
end

y = nan(length(v),nshuffles);
u = nanunique(c);
for j=1:nshuffles
    for i=1:length(u)
        inds = find(c==u(i));
        y(inds,j) = v(shuffle(inds));
    end
end