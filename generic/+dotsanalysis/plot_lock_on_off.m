function [p,Ton,Son,Toff,Soff,filter,lineas,out] = plot_lock_on_off(s,sm,conditions,filter,doPlot,...
    legend_str, varargin)
% function [p,Ton,Son,Toff,Soff,filter] = plot_lock_on_off(s,sm,conditions,filter,doPlot)

cutoff = 0.7;
independent_cutoff = false;
p_handle = [];
colortype = 0;
for i=1:2:length(varargin)
    if isequal(varargin{i},'cutoff')
        cutoff = varargin{i+1};
    elseif isequal(varargin{i},'independent_cutoff')
        independent_cutoff = varargin{i+1};
    elseif isequal(varargin{i},'p_handle')
        p_handle = varargin{i+1};
    elseif isequal(varargin{i},'color')
        colortype = varargin{i+1};
    elseif isequal(varargin{i},'colors')
        colores = varargin{i+1};
    end
end

if nargin<3 || isempty(filter)
    filter = ones(size(conditions,1),1)==1;
end

if nargin<5 || isempty(doPlot)
    doPlot = 1;
end
if nargin<6
    legend_str = [];
end

% if nargin<7
%     independent_cutoff = false;
% end

if doPlot == 2
    showErrorBars = true;
else
    showErrorBars = false;
end

lineas = [];
out = [];

Son  = s.g_on;
Soff = s.g_off;
Ton  = s.t_on;
Toff = s.t_off;

ncond = size(conditions,2);
if ~exist('colores','var')
    if ncond>1
        colores = rainbow_colors(ncond,'colortype',colortype);
        % colores = movshon_colors(ncond);
    else
        colores = [0 0 1];
    end
end


% cutoff = 0.7;
%ignore from cutoff rows that are all NaN
mn = mean(all(isnan(s.g_on')),2);
cutoff = mn + (1-mn)*cutoff;

if (doPlot>0)
    if ~isempty(p_handle)
        p = p_handle;
    else
        p = publish_plot(1,2);
    end
else
    p = [];
end

for i = 1:ncond
    ind = conditions(:,i)==1 & filter==1;
    [~,son,erron ] = curva_media(Son, ones(size(filter)) ,ind,0);
    [~,soff,erroff] = curva_media(Soff, ones(size(filter)),ind,0);
    
    doBootstrap = false;
    if (doBootstrap)
        [~,erron] = bootstr(Son(ind,:));
        [~,erroff] = bootstr(Soff(ind,:));
    end
    
    if (independent_cutoff)
        nans_on  = nanmean(isnan(Son(ind,:)));
        nans_off = nanmean(isnan(Soff(ind,:)));
    else
        nans_on  = nanmean(isnan(Son));
        nans_off = nanmean(isnan(Soff));
    end
    
    tind_on  = nans_on<cutoff;
    tind_off = nans_off<cutoff;
    
    
    
    I1 = tind_on & ~isnan(son);
    out(i).on.t=Ton(I1);
    out(i).on.y=son(I1);
    out(i).on.e=erron(I1);
    
    I2 = tind_off & ~isnan(soff);
    out(i).off.t=Toff(I2);
    out(i).off.y=soff(I2);
    out(i).off.e=erroff(I2);
    
    if (doPlot>0)
        p.next();
    
        if (showErrorBars)
            [~,l] = niceBars(Ton(I1),son(I1),erron(I1),colores(i,:),0.5);
            hold all
            set(l,'LineWidth',1);
            %for output
            
        else
            plot(Ton(tind_on),son(:,tind_on)','color',colores(i,:));
            hold all
        end

        p.next();

        if (showErrorBars)
            [~,lineas(i)] = niceBars(Toff(I2),soff(I2),erroff(I2),colores(i,:),0.5);
            hold all
            set(lineas,'LineWidth',1)

            

        else
            lineas(i) = plot(Toff(tind_off),soff(tind_off)','color',colores(i,:));
            hold all
        end
    end
end

if (doPlot>0)
    same_ylim(p.h_ax);
    
    if isempty(legend_str)
        % hl = legend_n(1:ncond,'hline',lineas);
        hl = [];
    else
        hl = legend(lineas,legend_str);
    end
    if ~isempty(hl)
        set(hl,'location','best','interpreter','none')
    end
    
    p.format('FontSize',15)
    
    %     set(hl,'FontSize',8)
    
    set(p.h_ax(2),'YAxisLocation','right')
    
    %ha = brokenAxis([0.5 0.11 0.012 0.06],'vertical');
    %set(gcf,'CurrentAxes',p.h_ax(2));
    
    set(gcf,'Position',[457  399  872  333])
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