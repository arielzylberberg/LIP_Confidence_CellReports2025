function z = zscore_bygroup(v,group)
%function z = zscore_bygroup(v,group)

if ~isvector(group)
    [~,~,group] = unique(group,'rows'); % overwrite "group" !!
end

transpose = false;
if ~isvector(v) && length(group)~=size(v,1)
    v = v';
    transpose = true;
end

ugroup = unique(group);
z = nan(size(v));
for i=1:length(ugroup)
    inds= group==ugroup(i);
    if isvector(v)
        z(inds) = nanzscore(v(inds));
    else
        temp = v(inds,:);
        ztemp = nanzscore(temp(:));
        z(inds,:) = reshape(ztemp,sum(inds),size(z,2));
    end
end

if transpose
    z = z';
end