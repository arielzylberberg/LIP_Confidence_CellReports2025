function [x,Y,S] = choice_against_duration(c,scoh,duration,nbins,plotflag,filt)
% function [x,Y] = choice_against_duration(c,scoh,duration,nbins,plotflag,filt)

if nargin<6 || isempty(filt)
    filt = ones(size(c))==1;
end

c = double(c);
c = c(filt);
scoh = scoh(filt);
duration = duration(filt);

coh = abs(scoh);
unicoh = unique(coh);
[x_aux,idx_prctile] = index_prctile(duration,linspace(0,100,nbins));
x = .5*x_aux(1:end-1) + .5*x_aux(2:end);
Y = nan(length(x),length(unicoh));
S = nan(length(x),length(unicoh));
for i=1:length(unicoh)
    inds = coh == unicoh(i);
    [~,aa,ss] = curva_media(c,idx_prctile,inds,0);
    ind = ismember(unique(idx_prctile(inds)),unique(idx_prctile));
    Y(ind,i) = aa;
    S(ind,i) = ss;
end

if plotflag>0
    rplot(x,Y,'Marker','.')
end

end