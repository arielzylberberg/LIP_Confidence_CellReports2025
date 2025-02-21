function [p,out] = do_plot_aligned_bigSbigT(id_neurons_to_include, RT, t, H, coh, include)

% select here

if ndims(H)==2 % ignore first argument if H has dim 2
    h = H;
else
    if (islogical(id_neurons_to_include) && sum(id_neurons_to_include)>1) || ...
        (~islogical(id_neurons_to_include) && length(id_neurons_to_include)>1)
        h = squeeze(nanmean(H(id_neurons_to_include,:,:)));
    else
        h = squeeze(H(id_neurons_to_include,:,:));%?
    end
end

%% new
RT = RT(include);
h = h(:,include);
coh = coh(include);
include = ones(length(RT),1)==1;

%%

h_aux = motionenergy.remove_post_decision_samples(h', t, RT - 0.1);

E = RT;
times = t;
data = h;
win = [5,5];
data_s = data;
% tind = t<0.25;
tind = t<0.2;
% tind = t<0.0;
% data_s(:,tind,:) = nan;
data_s(tind,:) = nan;
[datamp,timeslocked] = eventlockedmatc(data_s,times,E,win);

dt_s = t(2)-t(1);
s.t_on = t;
s.g_on = 1/dt_s * h_aux;


s.t_off = timeslocked;
s.g_off = 1/dt_s * datamp';


doPlot = 1;
filter = include;
conditions = adummyvar(coh);


% remove if all zero
ind = sum(conditions(filter,:))==0;
conditions(:,ind) = [];


colores = NS_colors(size(conditions,2));

% [p,Ton,Son,Toff,Soff,filter,lineas,out] = dotsanalysis.plot_lock_on_off(s,[],conditions,filter,...
%     doPlot,[],'independent_cutoff',1,'colors',colores);



smooth_win_t = 0.05;
% [p,lineas,out] = dotsanalysis.plot_lock_on_off_bigSbigT(s,smooth_win_t,conditions,filter,doPlot,...
%     [],'independent_cutoff',1,'colors',colores,'min_num_bins',500);

[p,lineas,out] = dotsanalysis.plot_lock_on_off_bigSbigT(s,smooth_win_t,conditions,filter,doPlot,...
    [],'independent_cutoff',1,'colors',colores,'min_num_bins',500,'cutoff',0.5);

% convert the out

Ton = out(1).on.t;
Toff = out(1).off.t;
for i=1:length(out)
    Son(i,:) = out(i).on.y;
    Soff(i,:) = out(i).off.y;
end



% 

set(p.h_ax(1),'xlim',[-0.2,1]);
% set(p.h_ax(1),'xlim',[-0.2,2]);
set(p.h_ax(2),'xlim',[-0.7,0.0]);
same_xscale(p.h_ax);

p.current_ax(1);
ylabel('Firing rate (sp/s)');
xlabel({'Time from','motion onset [s]'});

p.current_ax(2);
xlabel({'Time from','response [s]'});

p.shrink(1:2,1,0.9,1);
p.displace_ax(1:2,0.1,2);
p.format('FontSize',15,'LineWidthPlot',1);

clear out
out.Ton = Ton;
out.Son = Son;
out.Toff = Toff;
out.Soff = Soff;

end