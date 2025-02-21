function [centers,rth,rthc_norm] = hist_per_coh2(rt,coh,lim,n)
% function [centers,rth,rthc_norm] = hist_per_coh2(rt,coh,lim,n)
% Returns histograms of rt split by coh

edges = linspace(lim(1),lim(2),n);
ucoh = unique(coh);
ncoh = length(ucoh);
rth = nan(length(edges),ncoh);
for i=1:ncoh
    inds = coh==ucoh(i);
    rth(:,i) = histc(rt(inds),edges);
end

rth = rth(1:end-1,:);
centers = (edges(1:end-1)+edges(2:end))/2;

% cumsum, normalized
rthc_norm = bsxfun(@times,cumsum(rth,1),1./sum(rth,1));

% N = hist(Y,M)