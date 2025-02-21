function [X,Y,E] = choice_against_duration_var(c,scoh,duration,nbins,plotflag,filt)
% function [X,Y,E] = choice_against_duration_var(c,scoh,duration,nbins,plotflag,filt)

if nargin<6 || isempty(filt)
    filt = ones(size(c))==1;
end

c        = c(filt);
scoh     = scoh(filt);
duration = duration(filt);

coh = abs(scoh);
unicoh = unique(coh);

X = nan(nbins,length(unicoh));
Y = nan(nbins,length(unicoh));

for i=1:length(unicoh)
    inds = coh == unicoh(i);
    [x_aux,idx_prctile] = index_prctile(duration(inds),linspace(0,100,nbins+1));
    x = .5*x_aux(1:end-1) + .5*x_aux(2:end);
    
    X(:,i) = x;
    [~,Y(:,i),E(:,i)] = curva_media(c(inds),idx_prctile,[],0);
end
if plotflag>0
    if plotflag==2
    hp = errorbar(X,Y,E,'.-');
    else
        hp = plot(X,Y,'.-');
    end
    if length(hp)>1
        colores = rainbow_colors(length(hp),'colorType',5);
    else
        colores = [0,0,1];
    end
    for i=1:length(hp)
        set(hp(i),'Color',colores(i,:));
    end
    hl = legend_n(unicoh);
    set(hl,'location','best');
    
end