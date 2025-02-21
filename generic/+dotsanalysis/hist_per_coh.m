function [rth,edges,rtc_norm] = hist_per_coh(rt,coh,lim,n)

edges = linspace(lim(1),lim(2),n);
ucoh = unique(coh);
ncoh = length(ucoh);
rth = nan(length(edges),ncoh);
for i=1:ncoh
    inds = coh==ucoh(i);
    rth(:,i) = histc(rt(inds),edges);
end

% cumsum, normalized
rtc_norm = bsxfun(@times,cumsum(rth,1),1./sum(rth,1));