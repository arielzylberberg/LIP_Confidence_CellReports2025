function [C, x, uni_group] = hist_per_group(vals, group, ncats)

if nargin<3 || isempty(ncats)
    ncats = 40;
end

if isvector(group)
    group = group(:);
end

epsilon = nanstd(vals)/1000;
limi = [nanmin(vals(:))-epsilon, nanmax(vals(:))+epsilon];

edges = linspace(limi(1),limi(2), ncats);
dedges = [edges(2)-edges(1)]/2;
x = edges(1:end-1) + dedges/2;

uni_group = unique(group,'rows');

C = zeros(length(edges), length(uni_group));
for i=1:length(uni_group)
    I = all(group==uni_group(i,:),2);
    C(:,i) = histc(vals(I), edges);
end

C = C(1:end-1,:);

% rplot(x,C)

end
