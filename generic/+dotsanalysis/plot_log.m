function [hmark,hline] = plot_log(coh,ydata,coh_model,ymodel,varargin)
% function plot_log(coh,ydata,coh_model,ymodel,varargin)

color = 'b';
ylimit = nan;
break_scaling = 1;
marker = '.';
markersize = 25;
lsty = '-';
mfc = nan;
mec = nan;
show_errors = 1;
for i=1:length(varargin)
    if isequal(varargin{i},'color')
        color = varargin{i+1};
    elseif isequal(lower(varargin{i}),'linestyle')
        lsty = varargin{i+1};
    elseif isequal(varargin{i},'ylim')
        ylimit = varargin{i+1};
    elseif isequal(varargin{i},'break_scaling')
        break_scaling = varargin{i+1};
    elseif isequal(varargin{i},'marker')
        marker = varargin{i+1};
    elseif isequal(varargin{i},'markersize')
        markersize = varargin{i+1};
    elseif isequal(lower(varargin{i}),'markerfacecolor')
        mfc = varargin{i+1};
    elseif isequal(lower(varargin{i}),'markeredgecolor')
        mec = varargin{i+1};
    elseif isequal(lower(varargin{i}),'show_errors')
        show_errors = varargin{i+1};    
    end
end

if isnan(mfc)
    mfc = color;
end
if isnan(mec)
    mec = color;
end
%%
if nargin<3 || isempty(coh_model)
    coh_model = coh;
    ymodel = ydata;
end

%% data
if ~isempty(coh)
    [ucoh_data,uydata,stderr] = fold_for_logplot(coh,ydata);
end

%% model
if ~isempty(coh_model)
    [ucoh_model,uymodel] = fold_for_logplot(coh_model, ymodel,'model');
end

%%
hold on
if ~isempty(ucoh_model)
    hline = plot(ucoh_model,uymodel,'color',color,'LineWidth',0.75,'LineStyle',lsty);
else
    hline = [];
end
if ~isempty(coh)
    if show_errors
        terrorbar(ucoh_data,uydata,stderr,'color',color,'linestyle','none','marker','none');
        hmark = plot(ucoh_data,uydata,'color',color,'linestyle','none','marker',marker,...
             'markersize',markersize,'markerfacecolor',mfc,'markeredgecolor',mec);
    else
        hmark = plot(ucoh_data,uydata,'color',color,'linestyle','none','marker',marker,...
            'markersize',markersize,'markerfacecolor',mfc,'markeredgecolor',mec);
    end
end
set(gca,'xscale','log','xlim',[0.01,1])
% set(gca,'xtick',[0.016,0.032,0.064,0.128,0.256,0.512])
% set(gca,'xticklabel',{'0','3.2','6.4','12.8 ','25.6',' 51.2'})

set(gca,'xtick',[0.016,0.032,0.128,0.512])
set(gca,'xticklabel',{'0','3.2','12.8 ',' 51.2'})

if ~isnan(ylimit)
    ylim(ylimit);
end

% break axis
if break_scaling~=0
    BreakXAxis(gca,0.024,'break_scaling',break_scaling);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% HELPER FUNCTIONS %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [u,y,sterr] = fold_for_logplot(coh,v,type)

if nargin<3
    type = 'data';
end

if nanmax(coh)>1
    coh = coh/100;
end

acoh = abs(coh);
first_positive = 0.032;
zero_val = first_positive/2;

if isequal(type,'data')
    [u,y,sterr] = curva_media(v,acoh,[],0);
    sterr0 = sterr(u==0);
else
    [u,y] = curva_media(v,acoh,[],0);
    sterr = nan;
end
I = u<first_positive*0.75;
y(I) = nan;

% add zero and horizontal line

x2 = zero_val;
if isequal(type,'model')
    x3 = 1.2*zero_val;
    x1 = exp(2*log(x2)-log(x3));
    u = [x1; x2; x3; u(:)];
    val = nanmean(v(acoh==0));
    y = [val; val; val; y(:)];
else
    u = [x2; u(:)];
    val = nanmean(v(acoh==0));
    y = [val; y(:)];
    sterr = [sterr0; sterr(:)];
    
end

end

