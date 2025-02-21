function plot_raster(spk,t,varargin)

ntr = size(spk,1);
hline = false;
color = repmat([0,0,0],ntr,1);
for i=1:length(varargin)
    if isequal(varargin{i},'hline')
        hline = varargin{i+1};
    elseif isequal(varargin{i},'color')
        color = varargin{i+1};
    end
end


delta = 1; 

tmin = min(t);
tmax = max(t);
hold on
for i=1:ntr
    ind = spk(i,:)>0;
    f_ind = find(ind, 1);
    
    if sum(ind)>2
%         plot([t(ind)',t(ind)'],[i-0.4 i+0.4],'color',color(i,:))
        plot([t(ind)',t(ind)'],[i-delta i+delta],'color','k')
    elseif sum(ind)<=2
        tt=t(ind);
        for j=1:length(tt)
%             plot([tt(j),tt(j)],[i-0.4 i+0.4],'color',color(i,:))
            plot([tt(j),tt(j)],[i-delta i+delta],'color','k')
        end
    end
    if hline
        plot([tmin,tmax],[i,i],'k')
    end
end
ylim([-0.5,ntr+0.5])

end
