function [X,Y] = choice_against_duration_moving(c,scoh,duration,npoints,plotflag,filt)
% function [X,Y] = choice_against_duration_moving(c,scoh,duration,nbins,plotflag,filt)

if nargin<6 || isempty(filt)
    filt = ones(size(c))==1;
end

c = double(c);
c        = c(filt);
scoh     = scoh(filt);
duration = duration(filt);

% sort by duration
[~,inds] = sort(duration);
c        = c(inds);
scoh     = scoh(inds);
duration = duration(inds);
coh      = abs(scoh);
unicoh   = unique(coh);


X = cell(length(unicoh),1);
Y = cell(length(unicoh),1);
for i=1:length(unicoh)
    inds = coh == unicoh(i);
    x = conv(duration(inds),ones(npoints,1),'valid')/npoints;
    y = conv(c(inds),ones(npoints,1),'valid')/npoints;
    X{i} = x;
    Y{i} = y;
end

if plotflag>0
    n = max(3,length(unicoh));
    colores = rainbow_colors(n);
    for i=1:length(unicoh)
        plot(X{i},Y{i},'Marker','none','color',colores(i,:))
        hold all
    end
end

end