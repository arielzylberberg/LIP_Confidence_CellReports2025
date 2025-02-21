function [p,hlevel] = confidence_map(t,y,pcorrect,PHI,Bup,Blo,pcorrect_flag)

if nargin<7 || isempty(pcorrect_flag)
    pcorrect_flag = true;
end

logOdds = log10(pcorrect./(1-pcorrect));
pcorrect(pcorrect<0.5) = nan;
inds = mean(isnan(pcorrect),2)==1;
pcorrect(inds,:) = [];
logOdds(inds,:) = [];
y(inds) = [];

p = publish_plot(1,1);
if pcorrect_flag
    imagesc(t,y,pcorrect,[0.5,1])
else
    maxi = prctile(log10(pcorrect(:)./(1-pcorrect(:))),99.999);
    imagesc(t,y,logOdds,[0 maxi])
end


hold all
for i=1:length(PHI)
    phi = PHI(i);
    val = [log10(phi/(1-phi)),log10(phi/(1-phi))];
    levelCurves = contourcs(t,y,logOdds,val);
    [~,ind] = sort([levelCurves.Length],'descend');
    hlevel(i,1) = plot(levelCurves(ind(1)).X,levelCurves(ind(1)).Y,'w','LineWidth',2);
    hlevel(i,2) = plot(levelCurves(ind(2)).X,levelCurves(ind(2)).Y,'w','LineWidth',2);
end

lim = max(abs([Bup(:);Blo(:)])) + 0.1;

xlabel('time (s)')
ylabel('decision variable (a.u.)')

axis xy

hold on
plot(t,Bup,'k','LineWidth',2)
plot(t,Blo,'k','LineWidth',2)
xli = xlim;
xli(1)=0.03;
xlim(xli)
hc = colorbar;

ylim([-lim,lim])

if pcorrect_flag
    set(hc,'ytick',[0.5,1])
    set(get(hc,'title'),'string','p correct','FontSize',15)
else
    set(get(hc,'title'),'string','logodds correct','FontSize',15)
end


ht = text(0.2,0,'\Phi','interpreter','tex','color','w');

p.format('FontSize',25);
colormap('hot');


