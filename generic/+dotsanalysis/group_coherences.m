function coh_grouped = group_coherences(coh,group)

if nargin<2 || isempty(group)
    group = [1 1 2 2 ...
              2 3 4 ...
              4 4 5 5];
end

ucoh = nanunique(coh);
coh_grouped = nan(size(coh));
ugroup = unique(group);
for i=1:length(ugroup)
    inds = ismember(coh,ucoh(group==ugroup(i)));
    m = mean(ucoh(group==ugroup(i))); % assumes equal weight
    coh_grouped(inds) = m;
end