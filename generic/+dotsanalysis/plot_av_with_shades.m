function [p,lineas,out] = plot_av_with_shades(T,S,conditions,filter,doPlot,...
    legend_str,independent_cutoff)
% function [p,lineas,out] = plot_av_with_shades(T,S,conditions,filter,doPlot,legend_str,independent_cutoff)

if nargin<4 || isempty(filter)
    filter = ones(size(conditions,1),1);
end

if nargin<5 || isempty(doPlot)
    doPlot = 1;
end
if nargin<6
    legend_str = [];
end

if nargin<7
    independent_cutoff = false;
end

if doPlot == 2
    showErrorBars = true;
else
    showErrorBars = false;
end

lineas = [];
out = [];

% exand conditions if its just one column. TEST
% if size(conditions,2)==1
%     uu = nanunique(conditions);
%     cond = nan(length(conditions),length(uu));
%     for i=1:length(uu)
%         cond(:,i)=conditions==uu(i);
%     end
%     conditions = cond;
% end

ncond = size(conditions,2);
if ncond>1
    colores = rainbow_colors(ncond);
else
    colores = [0 0 1];
end


cutoff = 0.7;
%ignore from cutoff rows that are all NaN
mn = mean(all(isnan(S')),2);
cutoff = mn + (1-mn)*cutoff;

if (doPlot>0)
    p = publish_plot(1,1);
else
    p = [];
end

for i = 1:ncond
    ind = conditions(:,i)==1 & filter==1;
    [~,son,erron ] = curva_media(S, ones(size(filter)) ,ind,0);
    
    doBootstrap = false;
    if (doBootstrap)
        [~,erron] = bootstr(S(ind,:));
    end
    
    if (independent_cutoff)
        nans_on  = nanmean(isnan(S(ind,:)));
    else
        nans_on  = nanmean(isnan(S));
    end
    
    tind_on  = nans_on<cutoff;
    
    I1 = tind_on & ~isnan(son);
    out(i).t=T(I1);
    out(i).y=son(I1);
    out(i).e=erron(I1);
    
    
    if (doPlot>0)
        p.next();
    
        if (showErrorBars)
            [~,lineas(i)] = niceBars(T(I1),son(I1),erron(I1),colores(i,:),0.5);
            hold all
            set(lineas(i),'LineWidth',1);
            %for output
            
        else
            lineas(i) = plot(T(tind_on),son(:,tind_on)','color',colores(i,:));
            hold all
        end

    end
end

if (doPlot>0)
%     same_ylim(p.h_ax);
    
    if isempty(legend_str)
        hl = legend_n(1:ncond,'hline',lineas);
        set(hl,'location','best','interpreter','none')
        
    elseif ~iscell(legend_str) && legend_str==0
        % nothings
    else
        hl = legend(lineas,legend_str);
        set(hl,'location','best','interpreter','none')
    end
    
    
    p.format('FontSize',15)
    
    %     set(hl,'FontSize',8)
    
%     set(p.h_ax(2),'YAxisLocation','right')
    
    %ha = brokenAxis([0.5 0.11 0.012 0.06],'vertical');
    %set(gcf,'CurrentAxes',p.h_ax(2));
    
%     set(gcf,'Position',[457  399  872  333])
    set(gcf,'Position',[497  183  537  376])
end


function [mboot,sboot] = bootstr(y)
[n,nt] = size(y);
nboot = 100;
Sboot = zeros(nboot,nt);
for k=1:nboot %boot samples
    inds = ceil(rand(n,1)*n);
    Sboot(k,:) = nanmean(y(inds,:));
end
mboot = nanmean(Sboot);
sboot = nanstd(Sboot);