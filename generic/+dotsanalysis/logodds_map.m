function logodds_map(t,y,pcorrect,phi,hfig)

if nargin<5 || isempty(hfig)
    figure
else
    figure(hfig)
end

pcorrect = clip(pcorrect,.0001,1-.0001);
maxi = prctile(log10(pcorrect(:)./(1-pcorrect(:))),99.999);
logOdds = log10(pcorrect./(1-pcorrect));
imagesc(t,y,logOdds,[0 maxi])
hold all
val = [log10(phi/(1-phi)),log10(phi/(1-phi))];
levelCurves = contourcs(t,y,logOdds,val);
[~,ind] = sort([levelCurves.Length],'descend');
plot(levelCurves(ind(1)).X,levelCurves(ind(1)).Y,'k','LineWidth',2)
plot(levelCurves(ind(2)).X,levelCurves(ind(2)).Y,'k','LineWidth',2)
xlim([0 .8])
colorbar